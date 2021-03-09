#ifndef LINKSET_H
#define LINKSET_H

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <filesystem>
#include <fstream>  // for writing statistics to file
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <unordered_set>
#include <vector>

#include "mabl3/ProgressBar.h"
#include "mabl3/Timestep.h"
#include "nlohmann/json.hpp"
#include "prettyprint.hpp"
#include "tsl/hopscotch_map.h"
#include "tsl/hopscotch_set.h"
#include "Cube.h"
#include "DiagonalMatchesFilter.h"
#include "IdentifierMapping.h"
#include "Link.h"
#include "MemoryMonitor.h"
#include "SeedMap.h"

using namespace mabl3;



//! Helper function that computes the 'id'-th element of the cartesian product of a list of vectors
/*! \param id ID of the element to compute
 * \param c Vector of vectors, the inner vectors are the factors for the cartesian product
 *   (i.e. all possible combinations involving all inner vectors)
 *
 * \details For example, if the input is id=0 and c=[[1,2,3],[4,5,6],[7,8,9]], the result would
 *     be a vector [1,4,7]. For id=1, result is [2,4,7], id=3, result=[1,5,7] and so on.
 * */
template <typename T>
auto cartesianProductByID(size_t id,
                          std::vector<std::vector<T>> const & c) {
    std::vector<T> p;
    size_t outerProduct = 1;
    for (auto&& inner : c) {
        auto clock = std::floor(static_cast<double>(id) / static_cast<double>(outerProduct));
        auto innerID = static_cast<size_t>(clock) % inner.size();
        auto const & elem = *(std::next(inner.begin(), innerID));
        p.emplace_back(elem);
        outerProduct *= inner.size();
    }
    return p;
}



//! Class that can create and store Links of type LinkType
template<typename LinkType, typename LinkTypeHash, typename LinkTypeEqual>
class Linkset {
public:
    // ignore span: different masks may induce the same link with different spans, safe to ignore that
    using LinksetType = tsl::hopscotch_map<LinkType, size_t,
                                           LinkTypeHash, LinkTypeEqual>;
    //! Constructor (1)
    /*! \param config Shared ptr to global configuration
     * \param indetifierMapping Instance of IdentifierMapping that already
     * knows about all IDs
     * \param parallel Perform certain tasks in parallel if true
     *
     * \details Creates an empty Linkset with a user-defined IdentifierMapping
     * member that must already know all names of input genomes and sequences. */
    explicit Linkset(std::shared_ptr<Configuration const> config,
                     std::shared_ptr<IdentifierMapping const> identifierMapping,
                     bool parallel = true)
        : config_{config},
          idMapping_{identifierMapping},
          linkset_{}, numDiscarded_{0},
          parallel_{parallel && config->nThreads() > 1}, rd_{} {}
    //! Add a Link into the Linkset or increase counter for this Link
    void addLink(LinkType link);
    //! Applies M4, works only in 2D Case, with more dimensions only first two occurrences of Links used, may lead to UB
    void applyDiagonalMatchesFilter() {
        auto filter = DiagonalMatchesFilter<LinkType>(config_);
        Timestep ts("Filtering YASS-like");
        std::vector<LinkType> sortedMatches;
        for (auto&& elem : linkset_) { sortedMatches.emplace_back(elem.first); }
        linkset_.clear();
        if (parallel_) {
            std::sort(std::execution::par_unseq, sortedMatches.begin(), sortedMatches.end());
        } else {
            std::sort(sortedMatches.begin(), sortedMatches.end());
        }
        auto filteredMatches = filter.applyDiagonalMatchesFilter(sortedMatches);
        sortedMatches.clear();
        for (auto&& link : filteredMatches) { linkset_.insert({link, 0}); }
        filteredMatches.clear();
        std::cout << "[INFO] -- SeedMapSpaced::reportMatches -- Skipped " << filter.skippedNotInGenome1And2() << " matches that not included genome1 and 2" << std::endl;
        ts.endAndPrint();
    }
    //! Delete Link s in linkset
    void clear() {
        linkset_.clear();
        numDiscarded_ = 0;
    }
    //! Getter for member \c config_
    auto config() const { return config_; }
    //! Count number of valid links that can be created from this vector of occurrences
    size_t countLinks(std::vector<KmerOccurrence> const & occurrences) const {
        size_t count = 0;
        auto processingFunction = [this, &count](std::vector<tsl::hopscotch_map<size_t, // seqID
                                                                                tsl::hopscotch_set<KmerOccurrence,
                                                                                                   KmerOccurrencePositionHash,
                                                                                                   KmerOccurrencePositionEqual>>> const & occurrenceMap,
                                               size_t nPossible) {
            (void) occurrenceMap;
            if (nPossible > config_->matchLimit()) {
                count += config_->matchLimit();
            } else {
                count += nPossible;
            }
        };
        // count valid links
        processOccurrences(occurrences, processingFunction, config_->hasse());
        return count;
    }
    //! Create a Link in the Linkset from a vector of occurrences
    void createLinks(std::vector<KmerOccurrence> const & occurrences, size_t span);
    //! Create all Link s from a SeedMap
    template<typename TwoBitSeedDataType>
    void createLinks(SeedMap<TwoBitSeedDataType> const & seedMap, bool silent = false) {
        // parallel does slow things down, so use single thread
        ProgressBar pb(seedMap.size(), silent || config_->verbose() < 2);
        for (auto&& elem : seedMap.seedMap()) {
            for (size_t maskID = 0; maskID < config_->seedSetSize(); ++maskID) {
                auto& occurrenceVector = elem.second.at(maskID);
                createLinks(occurrenceVector, config_->maskCollection()->span(maskID));
            }
            ++pb;
        }
        pb.finish();
    }
    //! Create Link s from a single reference sequence vs. the other genomes
    template<typename TwoBitSeedDataType>
    void createLinks(SeedMap<TwoBitSeedDataType> const & seedMap, size_t sid) {
        if (seedMap.referenceSeedMap().referenceSeedMap().find(sid) == seedMap.referenceSeedMap().referenceSeedMap().end()) {
            throw std::runtime_error("[ERROR] -- Linkset::createLinks -- sid not found in referenceSeedMap");
        }
        if (idMapping_->sequenceIDToTuple().at(sid).gid != 0) {
            throw std::runtime_error("[ERROR] -- Linkset::createLinks -- sid not from reference genome");
        }
        for (auto&& seed : seedMap.referenceSeedMap().referenceSeedMap().at(sid)) {
            for (size_t maskID = 0; maskID < config_->seedSetSize(); ++maskID) {
                std::vector<KmerOccurrence> occurrenceVector{};
                for (auto&& occ : seedMap.seedMap().at(seed).at(maskID)) {
                    if (occ.genome() > 0 || occ.sequence() == sid) { occurrenceVector.emplace_back(occ); }
                }
                createLinks(occurrenceVector, config_->maskCollection()->span(maskID));
            }
        }
    }
    //! When metagraph as input, run this to create Links
/*    template<typename TwoBitSeedDataType>
    void createLinks(SeedKmerMap<TwoBitSeedDataType> const & seedKmerMap, bool silent = false) {
        (void) silent; // dummy for matching function signatures for SeedMap and SeedKmerMap
        auto processSeed = [this,
                            &seedKmerMap](TwoBitKmer<TwoBitSeedDataType> const & seed,
                                          Linkset<LinkType, LinkTypeHash, LinkTypeEqual> & linkset) {
            auto occurrenceVector = seedKmerMap.at(seed);
            for (size_t maskInd = 0; maskInd < config_->seedSetSize(); ++maskInd) {
                auto& occurrences = occurrenceVector.at(maskInd);
                linkset.createLinks(occurrences, config_->maskCollection()->span(maskInd));
            }
        };
        std::mutex mutex{};
        auto callback = [this,
                         &mutex,
                         &seedKmerMap,
                         &processSeed](typename SeedKmerMap<TwoBitSeedDataType>::SeedKmerMapType::const_iterator it,
                                       typename SeedKmerMap<TwoBitSeedDataType>::SeedKmerMapType::const_iterator end) {
           std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
           Linkset<LinkType, LinkTypeHash, LinkTypeEqual> linksetLocal{config_, idMapping_};
           for (; it != end; ++it) { processSeed(it->first, linksetLocal); }
           lock.lock();
           this->merge(linksetLocal);
           lock.unlock();
        };

        if (parallel_) {
            seedKmerMap.parallelIteration(callback);
        } else {
            for (auto&& elem : seedKmerMap.map()) { processSeed(elem.first, *this); }
        }
    }
    //! Create Link s from a single reference sequence vs. the other genomes (metagraph input)
    template<typename TwoBitSeedDataType>
    void createLinks(SeedKmerMap<TwoBitSeedDataType> const & seedKmerMap, size_t sid) {
        if (seedKmerMap.referenceSeedMap().referenceSeedMap().find(sid) == seedKmerMap.referenceSeedMap().referenceSeedMap().end()) {
            throw std::runtime_error("[ERROR] -- Linkset::createLinks -- sid not found in referenceSeedMap");
        }
        if (idMapping_->sequenceIDToTuple().at(sid).gid != 0) {
            throw std::runtime_error("[ERROR] -- Linkset::createLinks -- sid not from reference genome");
        }
        for (auto&& seed : seedKmerMap.referenceSeedMap().referenceSeedMap().at(sid)) {
            auto occurrenceVector = seedKmerMap.at(seed);
            for (size_t maskInd = 0; maskInd < config_->seedSetSize(); ++maskInd) {
                auto& occurrences = occurrenceVector.at(maskInd);
                createLinks(occurrences, config_->maskCollection()->span(maskInd));
            }
        }
    }
    //! Create Link s only in a predefined set of Cube s
    template<typename TwoBitSeedDataType>
    void createLinks(SeedKmerMap<TwoBitSeedDataType> const & seedKmerMap,
                     tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> const & relevantCubes,
                     bool silent = false) {
        (void) silent; // dummy for matching function signatures for SeedMap and SeedKmerMap
        auto processSeed = [this,
                            &seedKmerMap,
                            &relevantCubes](TwoBitKmer<TwoBitSeedDataType> const & seed,
                                            Linkset<LinkType, LinkTypeHash, LinkTypeEqual> & linkset) {
            auto occurrenceVector = seedKmerMap.at(seed);
            for (size_t maskInd = 0; maskInd < config_->seedSetSize(); ++maskInd) {
                auto& occurrences = occurrenceVector.at(maskInd);
                auto span = config_->maskCollection()->span(maskInd);
                linkset.createRelevantLinks(occurrences, span, relevantCubes);
            }
        };
        std::mutex mutex{};
        auto callback = [this,
                         &mutex,
                         &seedKmerMap,
                         &processSeed](typename SeedKmerMap<TwoBitSeedDataType>::SeedKmerMapType::const_iterator it,
                                       typename SeedKmerMap<TwoBitSeedDataType>::SeedKmerMapType::const_iterator end) {
           std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
           Linkset<LinkType, LinkTypeHash, LinkTypeEqual> linksetLocal{config_, idMapping_};
           for (; it != end; ++it) { processSeed(it->first, linksetLocal); }
           lock.lock();
           this->merge(linksetLocal);
           lock.unlock();
        };

        if (parallel_) {
            seedKmerMap.parallelIteration(callback);
        } else {
            for (auto&& elem : seedKmerMap.map()) { processSeed(elem.first, *this); }
        }
    }*/
    //! Create Link s only in a predefined set of Cube s
    template<typename TwoBitSeedDataType>
    void createLinks(SeedMap<TwoBitSeedDataType> const & seedMap,
                     tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> const & relevantCubes,
                     bool silent = false) {
        ProgressBar pb(seedMap.size(), silent || config_->verbose() < 2);
        for (auto&& elem : seedMap.seedMap()) {
            for (size_t maskID = 0; maskID < config_->seedSetSize(); ++maskID) {
                auto span = config_->maskCollection()->span(maskID);
                auto& occurrenceVector = elem.second.at(maskID);
                if (occurrenceVector.size()) { createRelevantLinks(occurrenceVector, span, relevantCubes); }
            }
            ++pb;
        }
        pb.finish();
    }
    //! Create a Link in the Linkset from a vector of occurrences
    void createRelevantLinks(std::vector<KmerOccurrence> const & occurrences, size_t span,
                             tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> const & relevantCubes);
    //! Group Links that overlap on the same diagonal
    void groupOverlappingLinks() {
        std::vector<LinkType> sortedLinks;
        for (auto&& elem : linkset_) { sortedLinks.emplace_back(elem.first); }
        linkset_.clear();
        if (parallel_) {
            std::sort(std::execution::par_unseq, sortedLinks.begin(), sortedLinks.end());
        } else {
            std::sort(sortedLinks.begin(), sortedLinks.end());
        }
        groupOverlappingSeeds(sortedLinks);
        for (auto&& link : sortedLinks) { linkset_.insert({link, 0}); }
    }
    //! Getter for the member variable \c idMapping_
    auto const & idMapping() const { return idMapping_; }
    //! Return number of times that \c link was created
    size_t linkCount(LinkType const & link) const {
        return (linkset_.find(link) != linkset_.end())
                ? linkset_.at(link)
                : 0;
    }
    //! Getter for the member variable \c linkset_
    auto const & linkset() const { return linkset_; }
    //! Merge \c rhs into this Linkset
    void merge(Linkset const & rhs) {
        numDiscarded_ += rhs.numDiscarded_;
        linkset_.insert(rhs.linkset_.begin(), rhs.linkset_.end());
    }
    //! Return number of discarded k-mers during construction of this Linkset
    size_t numDiscardedKmers() const { return numDiscarded_; }
    //! Return number of genomes in input data
    size_t numGenomes() const { return idMapping_->numGenomes(); }
    Linkset<LinkType, LinkTypeHash, LinkTypeEqual>& operator=(Linkset<LinkType, LinkTypeHash, LinkTypeEqual> && other) {
        config_ = other.config_;
        numDiscarded_ = other.numDiscarded_;
        idMapping_ = other.idMapping_;
        linkset_ = std::move(other.linkset_);
        parallel_ = other.parallel_;
        return *this;
    }
    //! Implements operator== for Linkset
    /*! Checks if \c linkset_ members are equal, i.e. the same Link s with the same counts
     * must be present, plus the \c idMapping_ members and the remaining members must be equal as well */
    bool operator==(Linkset<LinkType, LinkTypeHash, LinkTypeEqual> const & rhs) const {
        return linkset_ == rhs.linkset_
                && *idMapping_ == *(rhs.idMapping_)
                && config_->matchLimit() == rhs.config_->matchLimit()
                && config_->occurrencePerGenomeMax() == rhs.config_->occurrencePerGenomeMax()
                && config_->occurrencePerGenomeMin() == rhs.config_->occurrencePerGenomeMin();
    }
    //! Implements operator<< for Linkset objects for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, Linkset const & ls) {
        out << "ID Mapping:" << std::endl << *(ls.idMapping_) << std::endl;
        out << "Link set:" << std::endl;
        for (auto&& element : ls.linkset_) {
            out << "Count" << element.second << std::endl;
            out << element.first << std::endl;
        }
        return out;
    }
    std::ostream & printStatistics(std::ostream & out) {
        out << "Currently " << linkset_.size() << " seeds.\n";
        out << "Discarded " << numDiscarded_ << " seeds" << std::endl;
        return out;
    }
    bool processOccurrences(std::vector<KmerOccurrence> const & occurrenceVector,
                            std::function<void(std::vector<tsl::hopscotch_map<size_t, // seqID
                                                                              tsl::hopscotch_set<KmerOccurrence,
                                                                                                 KmerOccurrencePositionHash,
                                                                                                 KmerOccurrencePositionEqual>>> const &,
                                                size_t)> processingFunction,
                            bool hasse) const;
    //! Return number of Link s in this Linkset
    size_t size() const { return linkset_.size(); }

private:
    //! Configuration object
    std::shared_ptr<Configuration const> config_;
    //! Assigns IDs to genome and sequence strings in order of their appeareances
    /*! The reference genome always gets ID '0' */
    std::shared_ptr<IdentifierMapping const> idMapping_;
    //! Stores a mapping from Link s to their counts
    LinksetType linkset_;
    //! Count discarded or skipped k-mers
    size_t numDiscarded_;
    //! Perform certain tasks in parallel
    bool parallel_;
    //! Used to obtain seed for random link selection if there are too many possibilities
    std::random_device rd_;
};

#endif // LINKSET_H
