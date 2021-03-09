#ifndef SEEDKMERMAP_H
#define SEEDKMERMAP_H

#include <functional>
#include <memory>
#include <mutex>
#include <thread>

#include <tsl/hopscotch_map.h>
#include <mabl3/ProgressBar.h>
#include <MetagraphInterface.h>
#include "Configuration.h"
#include "ContainerChunks.h"
#include "ExtractSeeds.h"
#include "KmerOccurrence.h"
#include "ReferenceSeedMap.h"
#include "SeedMap.h"
#include "TwoBitKmer.h"



// Tag to mark this class as SeedKmerMap
struct SeedKmerMapTag{};



template<typename TwoBitSeedDataType>
class SeedKmerMapBase {
public:
    using ValueType = std::vector<std::vector<KmerOccurrence>>;
    using SeedKmerMapType = tsl::hopscotch_map<TwoBitKmer<TwoBitSeedDataType>,
                                               std::vector<MetagraphInterface::NodeID>,
                                               TwoBitKmerHash<TwoBitSeedDataType>>;
    using ChunkCallback = std::function<void(typename SeedKmerMapType::const_iterator it,
                                             typename SeedKmerMapType::const_iterator end)>;
    using ParallelCallback = std::function<void(TwoBitKmer<TwoBitSeedDataType> const &, ValueType const &)>;
    using IsSeedKmerMap = SeedKmerMapTag;

    //static struct PrefilterSeedMap{} prefilterSeedMap;
    SeedKmerMapBase(std::shared_ptr<Configuration const> config,
                    std::shared_ptr<IdentifierMapping const> idMap)
        : config_{config}, idMap_{idMap}, map_{},
          maskCollection_{config_->maskCollection()},
          referenceSeedMap_{*config}, stringhasher_{} {}
    SeedKmerMapBase(std::shared_ptr<Configuration const> config,
                    std::shared_ptr<IdentifierMapping const> idMap,
                    PrefilterSeedMapTag)
        : config_{config}, idMap_{idMap}, map_{},
          maskCollection_{config_->preMaskCollection()},
          referenceSeedMap_{*config}, stringhasher_{} {}
    //! Returns vector<vector<KmerOccurrence>> for seed, throws if seed not found
    ValueType at(TwoBitKmer<TwoBitSeedDataType> const & seed) const {
        ValueType value(config_->seedSetSize());
        tsl::hopscotch_map<MetagraphInterface::NodeID,
                           std::pair<std::string,
                                     std::vector<KmerOccurrence>>> occurrences{};
        for (auto&& nodeID : map_.at(seed)) {
            auto kmer = config_->metagraphInterface()->getKmer(nodeID);
            for (auto&& annotation : config_->metagraphInterface()->getAnnotation(nodeID)) {
                occurrences[nodeID].first = kmer;
                occurrences[nodeID].second.emplace_back(idMap_->queryGenomeIDConst(annotation.genome),
                                                        idMap_->querySequenceIDConst(annotation.sequence, annotation.genome),
                                                        annotation.bin_idx,
                                                        annotation.reverse_strand,
                                                        kmer);
            }
        }
        auto dummy = false;
        for (size_t m = 0; m < maskCollection_->masks().size(); ++m) {
            for (auto&& nodeID : map_.at(seed)) {
                auto& mask = maskCollection_->masks()[m];
                auto reconstructSeed = ExtractSeeds<TwoBitSeedDataType>::seedFromKmer(occurrences[nodeID].first, mask, dummy);
                if (reconstructSeed == seed) {
                    value[m].insert(value[m].end(), occurrences[nodeID].second.begin(), occurrences[nodeID].second.end());
                }
            }
        }
        return value;
    }
    void clear() noexcept { map_.clear(); }
    //! Remove seeds with only one occurrence per mask, DOES NOT CLEANUP REFERENCESEEDMAP!
    void cleanup() {
        std::vector<TwoBitKmer<TwoBitSeedDataType>> discard;
        std::mutex callbackMutex;
        auto callback = [this, &discard, &callbackMutex](typename SeedKmerMapType::const_iterator it,
                                                         typename SeedKmerMapType::const_iterator end) {
            std::unique_lock<std::mutex> lock(callbackMutex, std::defer_lock);
            std::vector<TwoBitKmer<TwoBitSeedDataType>> discardLocal;
            for (; it != end; ++it) {
                if (it->second.size() == 1
                        && config_->metagraphInterface()->graph()->get_labels(it->second[0]).size() == 1) {
                    // only one kmer induced this seed and this kmer only has one occurrence
                    discardLocal.emplace_back(it->first);
                } else {
                    // if multiple kmers induced this, but only one per mask, check if also only one occ per mask and discard
                    auto singleOccs = true;
                    for (auto&& occv : at(it->first)) {
                        if (occv.size() > 1) {
                            singleOccs = false;
                        }
                    }
                    if (singleOccs) { discardLocal.emplace_back(it->first); }
                }
            }
            lock.lock();
            discard.insert(discard.end(), discardLocal.begin(), discardLocal.end());
            lock.unlock();
        };
        parallelIteration(callback);
        // delete unneccessary seeds
        for (auto&& seed : discard) { map_.erase(seed); }
    }
    //! (Fwd) Remove duplicates in referenceSeedMap_ in parallel
    void cleanupReferenceSeeds(bool parallel = true) { referenceSeedMap_.cleanupReferenceSeeds(parallel); }
    //! Getter for config
    auto const & config() const { return config_; }
    //! Check if map is empty
    bool empty() const { return size() == 0; }
    //! Run seed extraction from metagraph
    void extractSeeds() {
        auto callback = [this](TwoBitKmer<TwoBitSeedDataType> seed, MetagraphInterface::NodeID nodeID, size_t) {
            map_[seed].emplace_back(nodeID);
            // PROPABLY BIG SLOWDOWN, FIX THIS IN THE FUTURE
            if (!config_->allvsall()) {
                for (auto&& annot : config_->metagraphInterface()->getAnnotation(nodeID)) {
                    auto gid = idMap_->queryGenomeIDConst(annot.genome);
                    if (gid == 0) { referenceSeedMap_.addSeed(idMap_->querySequenceIDConst(annot.sequence, annot.genome), seed); }
                }
            }
        };
        ExtractSeeds<TwoBitSeedDataType> extractor{maskCollection_, idMap_, config_};
        extractor.extractFromMetagraph(callback);
    }
    //! Getter for member idMap_
    auto const & idMap() const { return idMap_; }
    //! Getter for member map_
    auto const & map() const { return map_; }
    bool operator==(SeedKmerMapBase const & rhs) const {
        return config_ == rhs.config_
                && idMap_ == rhs.idMap_
                && map_ == rhs.map_;
    }
    //! Divide \c map_ in equal chunks for the available therads and execute \c callback for each chunk
    void parallelIteration(ChunkCallback const & callback) const {
        // spawn threads
        std::vector<std::thread> threads;
        for (size_t i = 0; i < config_->nThreads(); ++i) {
            auto chunkIt = getContainerChunkBegin(map_, i, config_->nThreads());
            auto chunkEnd = getContainerChunkEnd(map_, i, config_->nThreads());
            if (chunkIt != chunkEnd) {
                // https://thispointer.com/c11-start-thread-by-member-function-with-arguments/
                threads.push_back(std::thread(callback, chunkIt, chunkEnd));
            }
        }
        // wait until all threads are finished
        std::for_each(threads.begin(), threads.end(), [](std::thread & t) { t.join(); });
    }
    //! Call \c callback for each seed and its associated value, in parallel
    void parallelIteration(ParallelCallback const & callback) const {
        // Lambda to create value and call callback for the chunks
        auto processChunk = [this, &callback](typename SeedKmerMapType::const_iterator it,
                                              typename SeedKmerMapType::const_iterator end) {
            for (; it != end; ++it) {
                callback(it->first, at(it->first));
            }
        };
        parallelIteration(processChunk);
    }
    //! Getter for member \c referenceSeedMap_
    auto const & referenceSeedMap() const { return referenceSeedMap_; }
    //! Fwd size() member of map_
    auto size() const { return map_.size(); }

protected:
    std::shared_ptr<Configuration const> config_;
    std::shared_ptr<IdentifierMapping const> idMap_;
    SeedKmerMapType map_;
    //! Explicit pointer to masks (either standard or prefilter masks)
    std::shared_ptr<SpacedSeedMaskCollection const> maskCollection_;
    //! Stores seeds occurring in each reference sequence
    ReferenceSeedMap<TwoBitSeedDataType> referenceSeedMap_;
    //! For thinning
    std::hash<std::string> stringhasher_;
};



//! Iterator class to make SeedKmerMap iterable
template<typename TwoBitSeedDataType>
class SeedKmerMapIterator {
public:
    using iterator_category = std::forward_iterator_tag;//std::input_iterator_tag;
    using value_type = std::pair<TwoBitKmer<TwoBitSeedDataType> const,
                                 typename SeedKmerMapBase<TwoBitSeedDataType>::ValueType const>;
    using difference_type = value_type;
    using ponter = value_type const *;
    using reference = value_type const &;

    //! default c'tor
    SeedKmerMapIterator() : currentIt_{}, currentValue_{}, seedMapRef_{} {}
    //! c'tor -- Create iterator to seed
    /*! \param it Iterator to desired seed
     * \param map Reference to underlying SeedKmerMap
     * \details Initialize all members, setting iterator to it */
    SeedKmerMapIterator(typename SeedKmerMapBase<TwoBitSeedDataType>::SeedKmerMapType::const_iterator it,
                        SeedKmerMapBase<TwoBitSeedDataType> const & map)
        : currentIt_{it}, currentValue_{}, seedMapRef_{map} {}

    // satisfy EqualityComparable
    bool operator==(SeedKmerMapIterator const & rhs) const {
        return seedMapRef_ == rhs.seedMapRef_
                && currentIt_ == rhs.currentIt_;
    }
    bool operator!=(SeedKmerMapIterator const & rhs) const {
        return !(*this == rhs);
    }
    // satisfy dereferencable
    value_type const & operator*() {
        currentValue_ = std::make_shared<value_type>(currentIt_->first,
                                                     seedMapRef_.at(currentIt_->first));
        return *currentValue_;
    }
    value_type const * operator->() {
        (void)**this; // make sure currentValue_ is correct
        return &(*currentValue_);
    }
    // satisfy incrementable
    SeedKmerMapIterator & operator++() {
        ++currentIt_;
        return *this;
    }
    SeedKmerMapIterator operator++(int) {
        auto tmp = *this;
        ++*this;
        return tmp;
    }
    // satisfy Swappable
    void swap(SeedKmerMapIterator & rhs) {
        auto & seedMapRefBuffer = rhs.seedMapRef_;
        auto currentItBuffer = rhs.currentIt_;

        rhs.seedMapRef_ = seedMapRef_;
        rhs.currentIt_ = currentIt_;

        seedMapRef_ = seedMapRefBuffer;
        currentIt_ = currentItBuffer;
    }
    friend void swap(SeedKmerMapIterator & lhs, SeedKmerMapIterator & rhs) {
        lhs.swap(rhs);
    }

private:
    //! Iterator to current seed
    typename SeedKmerMapBase<TwoBitSeedDataType>::SeedKmerMapType::const_iterator currentIt_;
    //! For fake dereferencing
    /*! ptr to allow construction without value, shared to keep implicit copy, move etc. */
    std::shared_ptr<value_type const> currentValue_;
    //! Reference to underlying SeedKmerMapBase
    SeedKmerMapBase<TwoBitSeedDataType> const & seedMapRef_;
};



template<typename TwoBitSeedDataType>
class SeedKmerMap : public SeedKmerMapBase<TwoBitSeedDataType> {
public:
    using iterator_type = SeedKmerMapIterator<TwoBitSeedDataType>;
    using const_iterator_type = iterator_type;
    using TwoBitSeedDataType_ = TwoBitSeedDataType;

    SeedKmerMap(std::shared_ptr<Configuration const> config,
                std::shared_ptr<IdentifierMapping const> idMap)
        : SeedKmerMapBase<TwoBitSeedDataType>(config, idMap) {}
    SeedKmerMap(std::shared_ptr<Configuration const> config,
                std::shared_ptr<IdentifierMapping const> idMap,
                PrefilterSeedMapTag)
        : SeedKmerMapBase<TwoBitSeedDataType>(config, idMap, PrefilterSeedMapTag{}) {}

    auto begin() const { return iterator_type(this->map_.begin(), *this); }
    auto cbegin() const { return begin(); }
    auto cend() const { return iterator_type(this->map_.end(), *this); }
    auto end() const { return cend(); }
    //! Returns an iterator to the correct seed if seed was found, iterator to end() otherwise
    auto find(TwoBitKmer<TwoBitSeedDataType> const & seed) const {
        return iterator_type(this->map_.find(seed), *this);
    }
    //! Enable calls like `for (auto&& elem : seedKmerMap.seedMap())`
    auto const & seedMap() const { return *this; }
};

#endif // SEEDKMERMAP_H
