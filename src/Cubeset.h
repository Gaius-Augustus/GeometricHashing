#ifndef CUBESET_H
#define CUBESET_H

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <fstream>  // for writing statistics to file
#include <iostream>
#include <iterator>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>
#include "mabl3/ProgressBar.h"
#include "mabl3/Timestep.h"
#include "nlohmann/json.hpp"
#include "prettyprint.hpp"
#include "tsl/hopscotch_map.h"
#include "tsl/hopscotch_set.h"
#include "Cube.h"
#include "Linkset.h"
#include "MemoryMonitor.h"
#include "ParallelizationUtils.h"

using namespace mabl3;



//! Computes hash value of a Link, LinkPtr or Cube w.r.t. sequence properties
struct SequenceCombinationHash {
    template <typename T>
    size_t operator()(T const & sc) const {
        std::hash<size_t> hashfun;
        size_t seed = 0;
        for (size_t i = 0; i < sc.dimensionality(); ++i) {
            customCombineHash(seed, hashfun(sc.genome(i)));
            customCombineHash(seed, hashfun(sc.sequence(i)));
            customCombineHash(seed, hashfun(sc.reverse(i)));
        }
        return seed;
    }
};

//! Check if a Link, LinkPtr or Cube are equal w.r.t. sequence properties
struct SequenceCombinationEqual {
    template <typename T1, typename T2>
    bool operator()(T1 const & lhs, T2 const & rhs) const {
        if (lhs.dimensionality() != rhs.dimensionality()) { return false; }
        for (size_t i = 0; i < lhs.dimensionality(); ++i) {
            if (!(lhs.genome(i) == rhs.genome(i)
                  && lhs.sequence(i) == rhs.sequence(i)
                  && lhs.reverse(i) == rhs.reverse(i))) {
                return false;
            }
        }
        return true;
    }
};



//! Represents the set of all Cube s
/*! Cubes are generated from a Linkset and a mapping from
 * Cubes to a set of Link s that is contained in the respective cube,
 * as well as a mapping from Tiledistance s to the cubes that share these
 * distances is stored. Finally, artificial cubes and the relation between
 * cubes (Hasse diagram) are computed and stored. */
class Cubeset {
public:
    using CubeMapType = tsl::hopscotch_map<std::shared_ptr<Cube const>,
                                           tsl::hopscotch_set<LinkPtr, LinkPtrHash>,
                                           CubePtrHash, CubePtrEqual>;
    //using LinkCountType = tsl::hopscotch_map<LinkPtr, size_t,
    //                                         SequenceCombinationHash, SequenceCombinationEqual>;
    using ScoreMapType = tsl::hopscotch_map<double,
                                            std::vector<std::shared_ptr<Cube const>>>;
    using SubcubeMapType = tsl::hopscotch_map<std::shared_ptr<Cube const>,
                                              tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual>,
                                              CubePtrHash, CubePtrEqual>;
    using TileMapType = tsl::hopscotch_map<Tiledistance, std::vector<std::shared_ptr<Cube const>>,
                                           TiledistanceHash>;

    //! Factory function that computes the score of each cube
    /*! Executed during Hasse computation but after all respective
     * predecessors are known */
    double computeCubeScore(std::shared_ptr<Cube const> cube) const;
    double computeCubeScoreOld(std::shared_ptr<Cube const> cube) const;

    //! Constructor
    /*! \param linkset Linkset from which the Cubeset is computed
     *
     * \details Iterates over all Links s and creates the respective Cube s,
     * stores the mappings from Cube to set of Link s and from Tiledistance
     * to set of Cube s */
    Cubeset(std::shared_ptr<Linkset<LinkPtr,
                                    LinkPtrHashIgnoreSpan,
                                    LinkPtrEqualIgnoreSpan> const> linkset,
            std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> sequenceLengths,
            bool parallel = true);
    //! Constructor (2)
    /*! Throws as Linkset<Link> is not supported */
    Cubeset(std::shared_ptr<Linkset<Link,
                                    LinkHashIgnoreSpan,
                                    LinkEqualIgnoreSpan> const> linkset,
            std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> sequenceLengths,
            bool parallel = true)
        : cubeMap_{}, linkFraction_{0}, linkset_{},
          parallel_{parallel && linkset->config()->nThreads() > 1},
          sequenceLengths_{sequenceLengths},
          scoreToCube_{}, subcubeMap_{}, tilemap_{} {
        throw std::runtime_error("[ERROR] -- Cubeset::Cubeset (2) -- Link not supported, use LinkPtr");
    }

    //! Shortcut for Configuration stored in Linkset
    std::shared_ptr<Configuration const> config() const { return linkset_->config(); }
    //! Getter for member \c cubeAreaCutoff_
    auto cubeLengthCutoff() const  { return config()->cubeLengthCutoff(); }
    //! Getter for member \c cubeMap_
    auto const & cubeMap() const { return cubeMap_; }
    //! Group links in \c cube
    void groupLinksInCube(std::shared_ptr<Cube const> const & cube) {
        std::vector<LinkPtr> sortedLinks(cubeMap_.at(cube).begin(), cubeMap_.at(cube).end());
        cubeMap_.at(cube).clear();
        std::sort(sortedLinks.begin(), sortedLinks.end());
        groupOverlappingSeeds(sortedLinks);
        cubeMap_.at(cube).insert(sortedLinks.begin(), sortedLinks.end());
    }
    //! Group links in all cubes
    void groupLinksInCubes() {
        for (auto&& elem : cubeMap_) {
            groupLinksInCube(elem.first);
        }
    }
    //! Getter for idMap of underlying linkset
    auto const & idMap() const { return linkset_->idMapping(); }
    //! Returns \c std::unordered_set of Link s in \c cube
    auto const & links(std::shared_ptr<Cube const> const & cube) const { return cubeMap_.at(cube); }
    //! Returns set of Link s in \c cube and all sub-cubes
    tsl::hopscotch_set<LinkPtr,
                       LinkPtrHash> linksIncludingSubcubes(std::shared_ptr<Cube const> const & cube) const;
    //! Getter for member \c eta_
    auto normalizationParameter() const { return config()->cubeScoreNormalizationParameter(); }
    //! Getter for member \c nsubtiles_
    auto nsubtiles() const { return config()->cubeScoreParameter(); }
    //! Getter for member \c numGenomes_
    auto numGenomes() const { return linkset_->numGenomes(); }
    //! Getter for member \c sequenceLengths_
    auto const & sequenceLengths() const { return sequenceLengths_; }
    //! Getter for member \c scoreToCube_
    auto const & scoreToCube() const { return scoreToCube_; }
    //! Getter for member \c subcubes_
    auto const & subcubeMap() const { return subcubeMap_; }
    //! Getter for member \c tilemap_
    auto const & tilemap() const { return tilemap_; }
    //! Provide const access to the underlying Linkset
    auto const & underlyingLinkset() const { return linkset_; }
    //! Implements operator<< for a Cubeset object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, Cubeset const & cs) {
        out << "Cube => Set of Links" << std::endl;
        for (auto&& element : cs.cubeMap_) {
            out << element.first << " => " << element.second << std::endl;
        }
        out << "Score => Cubes" << std::endl;
        for (auto&& element : cs.scoreToCube_) {
            out << element.first << " => " << element.second << std::endl;
        }
        return out;
    }

private:
    //! Calculate chunk ID from link
    size_t getChunkID(Link const & link) const {
        if (config()->cubeScoreParameterChunks() > 0) {
            size_t sum = 0;
            for (auto&& occ : link.occurrence()) { sum += occ.position(); }
            return static_cast<size_t>( std::floor(static_cast<double>(sum) / static_cast<double>(config()->cubeScoreParameterChunks())) );
        } else {
            return 0ull;
        }
    }
    //! Get all true subcubes from \c cube, i.e. lower-dim cubes that share all tiledistances with \c cube
    std::vector<std::shared_ptr<Cube const>> getSubcubes(std::shared_ptr<Cube const> cube) const {
        std::vector<std::shared_ptr<Cube const>> subcubes;
        if (cube->dimensionality() <= 2) { return subcubes; }
        // each cube has reference (always first td), this always needs to match
        for (auto&& candidate : tilemap_.at(cube->tiledistance(0))) {
            if (cube->hasSubcube(*candidate)) {
                subcubes.emplace_back(candidate);
            }
        }
        return subcubes;
    }
    //! Get all true supercubes from \c cube, i.e. higher-dim cubes that include all tiledistances of \c cube
    std::vector<std::shared_ptr<Cube const>> getSupercubes(std::shared_ptr<Cube const> cube) const {
        std::vector<std::shared_ptr<Cube const>> supercubes;
        if (cube->dimensionality() == linkset_->idMapping()->numGenomes()) { return supercubes; }
        // each cube has reference (always first td), this always needs to match
        for (auto&& candidate : tilemap_.at(cube->tiledistance(0))) {
            if (candidate->hasSubcube(*cube)) {
                supercubes.emplace_back(candidate);
            }
        }
        return supercubes;
    }
    //! Stores the mapping from each Cube to a tuple of set of links
    /*! Implicitly stores the set of all Cube s as keys */
    CubeMapType cubeMap_;
    // ! Count links per sequence combination, abuse Link to define a unique sequence combination
    //LinkCountType linkCount_;
    //! Stores (number of links)/(number of possible links) for scoring
    long double linkFraction_;
    //! Linkset from which this Cubeset was created
    std::shared_ptr<Linkset<LinkPtr,
                            LinkPtrHashIgnoreSpan,
                            LinkPtrEqualIgnoreSpan> const> linkset_;
    //! Execute certain tasks in parallel
    bool parallel_;
    //! Map of sequence ID => sequence length
    std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> sequenceLengths_;
    //! Stores a score to Cube (s) mapping for all Cube s
    ScoreMapType scoreToCube_;
    //! Stores all true subcubes to a cube
    SubcubeMapType subcubeMap_;
    //! Map each present Tiledistance to all cubes that share this Tiledistance
    TileMapType tilemap_;
};



//! Store sids and cubes belonging to a sequence cluster
struct SequenceCluster {
    SequenceCluster() : sids{}, cubes{} {}
    tsl::hopscotch_set<size_t> sids;
    tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> cubes;
    friend std::ostream & operator<<(std::ostream & out, SequenceCluster const & c) {
        out << "\tSequence IDs in cluster: " << c.sids << std::endl;
        out << "\tCubes in cluster: " << c.cubes;
        return out;
    }
};



//! Finding relevant cubes with high weight, not storing links, only counting
class PrefilterCubeset {
public:
    // mapping cube to linkcount
    using CubeMapType = tsl::hopscotch_map<std::shared_ptr<Cube const>, size_t,
                                           CubePtrHash, CubePtrEqual>;
    using CubeSetType = tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual>;

    template <typename SeedMapType>
    PrefilterCubeset(SeedMapType const & seedMap,
                     std::shared_ptr<IdentifierMapping const> idMap,
                     std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> sequenceLengths,
                     std::shared_ptr<Configuration const> config,
                     bool silent = false)
        : config_{config}, cubeMap_{}, idMap_{idMap},
          nLinksTotal_{0}, relevantCubeSet_{},
          sequenceLengths_{sequenceLengths} {
        // count total links and fill cube count map
        if (!silent) { std::cout << "[INFO] -- counting links in cubes" << std::endl; }
        Timestep tsCount{"counting links in cubes", silent};
        ProgressBar pb(seedMap.size(), config_->verbose() < 2 || silent);
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset(config_, idMap_);
        // wrapper
        auto processingFunction = [this](
                std::vector<tsl::hopscotch_map<size_t,
                                               tsl::hopscotch_set<KmerOccurrence,
                                                                  KmerOccurrencePositionHash,
                                                                  KmerOccurrencePositionEqual>>> const & occurrenceMap,
                size_t nPossible) {
            countLinks(occurrenceMap, nPossible);
        };
        for (auto&& elem : seedMap.seedMap()) {
            for (auto&& occurrences : elem.second) {
                linkset.processOccurrences(occurrences, processingFunction, config_->preHasse());
            }
            ++pb;
        }
        tsCount.endAndPrint();
        // score and save relevant cubes
        if (!silent) { std::cout << "[INFO] -- scoring cubes via link count" << std::endl; }
        Timestep tsScore{"scoring cubes", silent};
        ProgressBar pbScore(cubeMap_.size(), config_->verbose() < 2 || silent);
        for (auto&& elem : cubeMap_) {
            if (elem.second >= config_->preLinkThreshold()) {
                // if preAddNeighbouringCubes, create neighbours to high-enough scoring cubes and add them as well
                if (config_->preAddNeighbouringCubes()) {
                    std::vector<std::vector<long long>> shifts;
                    size_t nneighbours = 1;
                    for (size_t i = 1; i < elem.first->dimensionality(); ++i) {
                        shifts.emplace_back(std::vector<long long>{-1,0,1});
                        nneighbours *= 3;
                    }
                    for (size_t i = 0; i < nneighbours; ++i) {
                        auto tdvector = elem.first->tiledistance();
                        auto dvector = cartesianProductByID(i, shifts);
                        std::vector<Tiledistance> newTdVector{tdvector.at(0)};
                        for (size_t j = 0; j < dvector.size(); ++j) {
                            Tiledistance newTd{tdvector.at(j+1).genome(),
                                               tdvector.at(j+1).sequence(),
                                               tdvector.at(j+1).distance() + dvector.at(j),
                                               tdvector.at(j+1).reverse()};
                            newTdVector.emplace_back(newTd);
                        }
                        if (newTdVector.size() != tdvector.size()) { throw std::runtime_error("[ERROR] -- Neighbour creation failed"); }
                        if (newTdVector.at(0) != tdvector.at(0)) { throw std::runtime_error("[ERROR] -- Neighbour creation failed"); }
                        auto cube = std::make_shared<Cube>(newTdVector);
                        relevantCubeSet_.insert(cube);
                        // if Hasse in 2nd run, also create genome1-genome2 cubes for possible 2D links
                        if (cube->dimensionality() > 2 && config_->hasse()) {
                            auto lowerDimCube = std::make_shared<Cube>(std::vector<Tiledistance>{
                                                                           cube->tiledistance(0),
                                                                           cube->tiledistance(1)
                                                                       });
                            relevantCubeSet_.insert(lowerDimCube);
                        }
                    }
                } else {
                    // just add the high-enough scoring cube and possibly lower-dim parts
                    relevantCubeSet_.insert(elem.first);
                    // if Hasse in 2nd run, also create genome1-genome2 cubes for possible 2D links
                    if (elem.first->dimensionality() > 2 && config_->hasse()) {
                        auto lowerDimCube = std::make_shared<Cube>(std::vector<Tiledistance>{
                                                                       elem.first->tiledistance(0),
                                                                       elem.first->tiledistance(1)
                                                                   });
                        relevantCubeSet_.insert(lowerDimCube);
                    }
                }
            }
            ++pbScore;
        }
        tsScore.endAndPrint();
        if (!silent) { std::cout << "[INFO] -- found " << relevantCubeSet_.size() << " relevant cubes" << std::endl; }
    }

    auto const & idMap() const {return idMap_; }
    auto const & cubeMap() const { return cubeMap_; }
    auto const & relevantCubeSet() const { return relevantCubeSet_; }
    //! Create sequence clusters
    /*! Inside a cluster, there are relevant cubes including (a subset of) these sequences.
     *  There are no relevant cubes that include sequences from two or more different clusters */
    auto sequenceCluster() const {
        std::vector<std::shared_ptr<SequenceCluster>> sidCluster; // each set is a cluster of sids
        tsl::hopscotch_map<size_t, std::shared_ptr<SequenceCluster>> sidToCluster;
        for (auto&& cube : relevantCubeSet_) {
            bool belongsToCluster = false;
            std::shared_ptr<SequenceCluster> currentCluster;
            for (auto&& td : cube->tiledistance()) {
                if (sidToCluster.find(td.sequence()) != sidToCluster.end()) {
                    belongsToCluster = true;
                    currentCluster = sidToCluster.at(td.sequence());
                    break;
                }
            }
            if (!belongsToCluster) {
                // create new cluster
                currentCluster = std::make_shared<SequenceCluster>();
                sidCluster.emplace_back(currentCluster);
            }
            // make sure the cube and all sequences in cube get added and point to the cluster
            currentCluster->cubes.insert(cube);
            for (auto&& td : cube->tiledistance()) {
                currentCluster->sids.insert(td.sequence());
                sidToCluster[td.sequence()] = currentCluster;
            }
        }
        return sidCluster;
    }
    auto sequenceTuplesToCubes() const {
        tsl::hopscotch_map<std::vector<size_t>,
                           tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual>,
                           VectorHash> tuples;
        for (auto&& cube : relevantCubeSet_) {
            std::vector<size_t> tuple;
            for (auto&& td : cube->tiledistance()) {
                tuple.emplace_back(td.sequence());
            }
            tuples[tuple].insert(cube);
        }
        return tuples;
    }

private:
    //! counts total number of observed links from an occurrenceMap (callback for linkset.processOccurrences)
    void countLinks(std::vector<tsl::hopscotch_map<size_t, tsl::hopscotch_set<KmerOccurrence,
                                                                              KmerOccurrencePositionHash,
                                                                              KmerOccurrencePositionEqual>>> const & occurrenceMap,
                                                   size_t nPossible) {
        // create all links (resp. cubes) and count links per cube
        nLinksTotal_ += nPossible;
        // flatten occurrenceMap
        std::vector<std::vector<KmerOccurrence>> allOccs{};
        for (size_t i = 0; i < idMap_->numGenomes(); ++i) {
            std::vector<KmerOccurrence> occs;
            for (auto&& elem : occurrenceMap.at(i)) {
                if (elem.second.size()) {
                    occs.insert(occs.end(), elem.second.begin(), elem.second.end());
                }
            }
            if (occs.size()) {
                allOccs.emplace_back(occs);
            }
        }
        // count cubes
        for (size_t id = 0; id < nPossible; ++id) {
            Link link{cartesianProductByID(id, allOccs), 1}; // span doesn't matter
            auto cube = std::make_shared<Cube>(link, config_->tileSize());
            if (cubeMap_.find(cube) == cubeMap_.end()) {
                cubeMap_[cube] = 0;
            }
            cubeMap_.at(cube) += 1;
        }
    }

    std::shared_ptr<Configuration const> config_;
    CubeMapType cubeMap_;
    std::shared_ptr<IdentifierMapping const> idMap_;
    size_t nLinksTotal_;
    CubeSetType relevantCubeSet_;
    std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> sequenceLengths_;
};

#endif // CUBESET_H
