#ifndef SEEDMAP_H
#define SEEDMAP_H

#include <algorithm>
#include <cstdlib>
#include <execution>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <type_traits>

#include "mabl3/JsonStream.h"
#include "mabl3/ProgressBar.h"
#include "tsl/hopscotch_map.h"
#include "tsl/hopscotch_set.h"
#include "ContainerChunks.h"
#include "ExtractSeeds.h"
#include "IdentifierMapping.h"
#include "KmerOccurrence.h"
#include "ParallelizationUtils.h"
#include "ReferenceSeedMap.h"
#include "SpacedSeedMask.h"
#include "TwoBitKmer.h"
 
using namespace mabl3;

struct PrefilterSeedMapTag{};

//! Base class for a (spaced) seed map
/*! Dictates the methods a productive child needs to implement in order to work,
 * implements most methods */
template <typename TwoBitSeedDataType>
class SeedMap {
public:
    using TwoBitSeedDataType_ = TwoBitSeedDataType;
    // https://tessil.github.io/hopscotch-map/classtsl_1_1hopscotch__map.html
    // template<class Key, class T, class Hash = std::hash<Key>, class KeyEqual = std::equal_to<Key>,
    //          class Allocator = std::allocator<std::pair<Key, T>>, unsigned int NeighborhoodSize = 62,
    //          bool StoreHash = false, class GrowthPolicy = tsl::hh::power_of_two_growth_policy<2>>
    using SeedMapType = tsl::hopscotch_map<TwoBitKmer<TwoBitSeedDataType>,
                                           std::vector<std::vector<KmerOccurrence>>,    // outer vector: one element for each mask
                                           TwoBitKmerHash<TwoBitSeedDataType>,
                                           // default template parameters, except the last (growth policy with smallest allowed growth factor of 1.1)
                                           std::equal_to<TwoBitKmer<TwoBitSeedDataType>>,
                                           std::allocator<std::pair<TwoBitKmer<TwoBitSeedDataType>, std::vector<std::vector<KmerOccurrence>>>>,
                                           62, false, tsl::hh::mod_growth_policy<std::ratio<11,10>>>;

    //! Sanity check idMap during construction
    void factory() {
        if ((!idMap_->genomeKnown(config_->genome1())) || (idMap_->queryGenomeIDConst(config_->genome1()) != 0)) {
            throw std::runtime_error("[ERROR] -- SeedMap -- '--genome1' must have ID '0'");
        }
        if ((!idMap_->genomeKnown(config_->genome2())) || (idMap_->queryGenomeIDConst(config_->genome2()) != 1)) {
            throw std::runtime_error("[ERROR] -- SeedMap -- '--genome2' must have ID '1'");
        }
    }
    //! c'tor (1)
    /*! \param config Shared main configuration object
     * \param idMap Shared identifier map
     * \details Initializes all members */
    SeedMap(std::shared_ptr<Configuration const> config,
            std::shared_ptr<IdentifierMapping const> idMap)
        : config_{config}, idMap_{idMap},
          maskCollection_{config_->maskCollection()},
          mutex_{}, rd_{}, referenceSeedMap_{*config},
          seedMap_{} { factory(); }
    //! c'tor (2)
    /*! \param config Shared main configuration object
     * \param idMap Shared identifier map
     * \details Initializes all members, using seeds for prefilter-step in GH */
    SeedMap(std::shared_ptr<Configuration const> config,
            std::shared_ptr<IdentifierMapping const> idMap,
            PrefilterSeedMapTag)
        : config_{config}, idMap_{idMap},
          maskCollection_{config_->preMaskCollection()},
          mutex_{}, rd_{}, referenceSeedMap_{*config},
          seedMap_{} { factory(); }
    //! Copy c'tor
    /*! Does not copy mutexes and rng state */
    SeedMap(SeedMap<TwoBitSeedDataType> const & other)
        : config_{other.config_}, idMap_{other.idMap_},
          maskCollection_{other.maskCollection_},
          mutex_{}, rd_{}, referenceSeedMap_{other.referenceSeedMap_},
          seedMap_{other.seedMap_} {}

    //! Add a new seed to the seed map
    /*! This includes filtering and calling cleanup if neccessary, it is important that the
     * \c seedMap_ member is in a valid state according to filter criteria after calling \c addSeed() */
    void addSeed(TwoBitKmer<TwoBitSeedDataType> const & seed,
                 KmerOccurrence const & occurrence,
                 size_t maskIndex) {
        if (seed.length() != maskCollection_->weight()) { throw std::runtime_error("[ERROR] -- SeedMap::addSeed() -- Wrong seed length"); }
        accessSeedInMap(seed).at(maskIndex).emplace_back(occurrence);
        if ((!config_->allvsall()) && occurrence.genome() == 0) { referenceSeedMap_.addSeed(occurrence.sequence(), seed); }
    }
    //! Getter for member config_
    auto config() const { return config_; }
    //! (Fwd) Remove duplicates in referenceSeedMap_ in parallel
    void cleanupReferenceSeeds(bool parallel = true) { referenceSeedMap_.cleanupReferenceSeeds(parallel); }
    //! Clear the seedMap member that stores the mapping from a seed to its occurrences
    void clear() { seedMap_.clear(); }
    //! Run seed extraction from input fastas
    void extractSeeds(std::shared_ptr<FastaCollection const> fastaCollection,
                      bool parallel = true) {
        auto wrapper = [this](TwoBitKmer<TwoBitSeedDataType> seed, KmerOccurrence occ, size_t maskInd) {
            addSeed(seed, occ, maskInd);
        };

        ExtractSeeds<TwoBitSeedDataType> extractor{maskCollection_, idMap_, config_};
        extractor.extractFromFastas(fastaCollection, wrapper, parallel);
        cleanupReferenceSeeds(parallel);
    }
    //! Getter for member \c genome0_
    auto const & genome0() const { return config_->genome1(); }
    //! Getter for member \c genome1_
    auto const & genome1() const { return config_->genome2(); }
    //! Getter for member \c idMap_
    auto idMap() const { return idMap_; }
    //! Returns a copy of \c idMap_
    auto idMapCopy() const { return *idMap_; }
    //! Forward getter for member \c spacedSeedMasks_->masks()
    auto const & masks() const { return maskCollection_->masks(); }
    //! Merge \c seedMap_ of \c localMap with this
    void merge(std::shared_ptr<SeedMap<TwoBitSeedDataType>> localMap) {
        if (maskCollection_ != localMap->maskCollection_) { throw std::runtime_error("[ERROR] -- SeedMap::merge() -- Merging from different mask collections"); }
        for (auto&& elem : localMap->seedMap_) {
            auto& seed = elem.first;
            for (size_t i = 0; i < maskCollection_->size(); ++i) {
                auto& localOccurrenceVector = elem.second.at(i);
                accessSeedInMap(seed).at(i).insert(seedMap_.at(seed).at(i).end(),
                                                   localOccurrenceVector.begin(), localOccurrenceVector.end());
            }
        }
        referenceSeedMap_.merge(localMap->referenceSeedMap());
    }
    //! Forward getter for nThreads() in config
    auto nThreads() const { return config_->nThreads(); }
    //! Forwards to method \c numGenomes of \c idMap_
    auto numGenomes() const { return idMap_->numGenomes(); }
    //! Forwards to method \c numSequences of \c idMap_
    auto numSequences() const { return idMap_->numSequences(); }
    //! Copy assignment
    SeedMap<TwoBitSeedDataType> & operator=(SeedMap<TwoBitSeedDataType> const & other) {
        config_ = other.config_;
        idMap_ = other.idMap_;
        maskCollection_ = other.maskCollection_;
        referenceSeedMap_ = other.referenceSeedMap_;
        seedMap_ = other.seedMap_;
        return *this;
    }
    //! Print statistics about the seed map creation process
    void printStatistics() const {
        auto numKmers = seedMap_.size();
        std::cout << "Number of unique k-mers created from the input files: " << numKmers << std::endl;
    }
    //! Getter for member \c referenceSeedMap_
    auto const & referenceSeedMap() const { return referenceSeedMap_; }
    //! Getter for member \c seedSetSize_
    auto seedSetSize() const { return maskCollection_->size(); }
    //! Getter for member \c seedMap_
    auto const & seedMap() const { return seedMap_; }
    //! Return number of seeds
    auto size() const { return seedMap_.size(); }
    //! Getter for member spacedSeedMasks_
    auto spacedSeedMasks() const { return maskCollection_; }
    //! Getter for member \c span_
    auto span() const { return maskCollection_->maxSpan(); }
    //! Getter for member \c weight_
    auto weight() const { return maskCollection_->weight(); }

protected:
    //! Returns reference to \c seedMap_[seed], creates correctly sized vectors if new seed
    auto& accessSeedInMap(TwoBitKmer<TwoBitSeedDataType> const & seed) {
        if (seedMap_.find(seed) == seedMap_.end()) {
            seedMap_.emplace(seed, maskCollection_->size());
        }
        return seedMap_[seed];
    }

    //! Main configuration
    std::shared_ptr<Configuration const> config_;
    //! Identifier Mapping
    std::shared_ptr<IdentifierMapping const> idMap_;
    //! Explicit pointer to masks (either standard or prefilter masks)
    std::shared_ptr<SpacedSeedMaskCollection const> maskCollection_;
    //! General lock
    std::mutex mutex_;
    //! Used to obtain seed for random link selection if there are too many possibilities
    std::random_device rd_;
    //! Stores seeds occurring in each reference sequence
    ReferenceSeedMap<TwoBitSeedDataType> referenceSeedMap_;
    //! seedMap
    SeedMapType seedMap_;
};

#endif // SEEDMAP_H
