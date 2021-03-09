#ifndef REFERENCESEEDMAP_H
#define REFERENCESEEDMAP_H

#include <algorithm>

#include "Configuration.h"
#include "ParallelizationUtils.h"
#include "TwoBitKmer.h"

template<typename TwoBitSeedDataType>
class ReferenceSeedMap {
public:
    // store seeds occurring in each reference (genome0) sequence
    using ReferenceSeedMapType = tsl::hopscotch_map<size_t, std::vector<TwoBitKmer<TwoBitSeedDataType>>>;

    //! c'tor
    ReferenceSeedMap(Configuration const & config) : nThreads_{config.nThreads()}, referenceSeedMap_{} {}

    //! Add a single seed to sid, no check if sid is actually valid!
    void addSeed(size_t sid, TwoBitKmer<TwoBitSeedDataType> const & seed) {
        referenceSeedMap_[sid].emplace_back(seed);
    }
    //! Remove duplicates in referenceSeedMap_ in parallel
    void cleanupReferenceSeeds(bool parallel = true) {
        auto callback = [](ReferenceSeedMapType & referenceSeedMap,
                           typename ReferenceSeedMapType::const_iterator it,
                           typename ReferenceSeedMapType::const_iterator end) {
            for (; it != end; ++it) {
                auto& seedVector = referenceSeedMap.at(it->first);
                std::sort(seedVector.begin(), seedVector.end());
                auto last = std::unique(seedVector.begin(), seedVector.end());
                seedVector.erase(last, seedVector.end());   // remove indeterminate elements resulting from std::unique
                seedVector.shrink_to_fit();
            }
        };
        size_t nThreads = (parallel) ? nThreads_ : 1;
        executeParallel(referenceSeedMap_, nThreads, callback, std::ref(referenceSeedMap_));
    }
    //! Forward clear method of referenceSeedMap_
    void clear() { referenceSeedMap_.clear(); }
    //! Add all seeds from \c rhs to this map
    void merge(ReferenceSeedMap const & rhs) {
        for (auto&& elem : rhs.referenceSeedMap_) {
            auto& seedVector = referenceSeedMap_[elem.first];
            seedVector.insert(seedVector.end(), elem.second.begin(), elem.second.end());
        }
    }
    //! Getter for \c referenceSeedMap_
    auto const & referenceSeedMap() const { return referenceSeedMap_; }

private:
    size_t nThreads_;
    ReferenceSeedMapType referenceSeedMap_;
};

#endif // REFERENCESEEDMAP_H
