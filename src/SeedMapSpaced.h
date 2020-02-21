#ifndef SEEDMAPSPACED_H
#define SEEDMAPSPACED_H

#include <array>
#include <algorithm>
#include <mutex>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "boost/dynamic_bitset.hpp"
#include "Configuration.h"
#include "CustomHashGeneral.h"
#include "DiagonalMatchesFilter.h"
#include "KmerOccurrence.h"
#include "SeedMapGeneral.h"
#include "SpacedSeedMask.h"
#include "Timestep.h"
#include "TwoBitKmer.h"

//! Implements SeedMapGeneral methods shared among the SpacedSeed classes
template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
class SeedMapSpaced : public SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType> {
public:
    using MatchesType = std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash>;

    //! c'tor -- unified SeedMap child c'tor
    /*! \param config Configuration object
     * \param idMap Global identifier mapping
     * \param masks Collection of \c seedSetSize distinct spaced seed masks of length \c span and \c weight care-positions
     *
     * \details Sets everything up so that seeds from ExtractSeeds can be fed into \c this */
    SeedMapSpaced(std::shared_ptr<Configuration const> config,
                  std::shared_ptr<IdentifierMapping> idMap,
                  std::shared_ptr<SpacedSeedMaskCollection const> masks)
        : SeedMapGeneral<TwoBitKmerDataType,
                         TwoBitSeedDataType>(masks,
                                             config->genome1(),
                                             config->genome2(),
                                             idMap,
                                             config->matchLimit(),
                                             config->matchLimitDiscardSeeds(),
                                             config->createAllMatches(),
                                             config->nThreads(),
                                             config->quiet()),
          config_{config}, matches_{} {}

    void appendMatches(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> sourceMap) {
        auto sourceCast = std::static_pointer_cast<SeedMapSpaced<TwoBitKmerDataType, TwoBitSeedDataType>>(sourceMap);
        appendMatches(sourceCast->matches());
    }
    void appendMatches(MatchesType const & sourceMatches) {
        matches_.insert(sourceMatches.begin(),
                        sourceMatches.end());
    }
    void applyDiagonalMatchesFilter(std::shared_ptr<std::array<size_t, 4>> skipped);
    void clearMatches() { matches_.clear(); }
    void createMatches() {
        this->createMatchesGeneral();
    }
    void insertMatch(KmerOccurrence const & occ1, KmerOccurrence const & occ2) {
        matches_.emplace(occ1, occ2);
    }
    std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> localCopy() const;
    //! Getter for member \c matches_
    auto const & matches() const { return matches_; }
    void merge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType> > localMap);
    //! Merge matches from a source to member \c matches_
    template <typename ContainerSource>
    void mergeMatches(ContainerSource const & matchSource) {
        this->mergeMatches_impl(matchSource, matches_);
    }
    size_t numMatches() const { return matches_.size(); }
    auto orderedMatches() const {
        std::set<KmerOccurrencePair> orderedMatches;
        for (auto&& match : matches_) { orderedMatches.emplace(match); }
        return std::move(orderedMatches);
    }
    void output(std::ostream & outstream) const;
    void quickMerge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType> > localMap);

protected:
    //! Configuration for filter
    std::shared_ptr<Configuration const> config_;
    //! Store matches
    MatchesType matches_;
};

#endif // SEEDMAPSPACED_H
