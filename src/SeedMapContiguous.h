#ifndef EXACTMATCHES_H
#define EXACTMATCHES_H

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

#include "json/json.hpp"
#include "prettyprint/prettyprint.hpp"
#include "Configuration.h"
#include "KmerOccurrence.h"
#include "SeedMapGeneral.h"
#include "SpacedSeedMask.h"
#include "TwoBitKmer.h"

//! Productive class to execute method M1
template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
class SeedMapContiguous : public SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType> {
public:
    using MatchesType = std::vector<KmerOccurrencePair>;

    //! c'tor (1)
    /*! \param span k-mer length
     * \param genome0 Name of first genome
     * \param genome1 Name of second genome
     * \param idMap Global identifier mapping
     * \param globalOccurrenceThreshold Discard k-mers that globally occur more often than this
     * \param perSequenceOccurrenceThreshold Discard k-mers that occur more often than this in any sequence
     * \param createAllMatches In Triangulation Filter setting, all matches are required, so set this; otherwise can be false to save ressources
     * \param nThreads Number of threads for parallelization
     * \param quiet If true, do not print progress bars
     *
     * \details Sets everything up so that seeds from ExtractSeeds can be fed into \c this */
    SeedMapContiguous(size_t span, std::string const & genome0, std::string const & genome1,
                 std::shared_ptr<IdentifierMapping> idMap,
                 size_t matchLimit, bool matchLimitDiscardSeeds, bool createAllMatches,
                 size_t nThreads, bool quiet = false)
        : SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>(span, genome0, genome1, idMap,
                                                                 matchLimit, matchLimitDiscardSeeds,
                                                                 createAllMatches,
                                                                 nThreads, quiet),
          matches_{} {}
    //! c'tor (2) -- unified SeedMap child c'tor
    /*! \param config Configuration object
     * \param idMap Global identifier mapping
     * \param masks Collection of \c seedSetSize distinct spaced seed masks of length \c span and \c weight care-positions, not use here!
     *
     * \details Sets everything up so that seeds from ExtractSeeds can be fed into \c this */
    SeedMapContiguous(std::shared_ptr<Configuration const> config,
                      std::shared_ptr<IdentifierMapping> idMap,
                      std::shared_ptr<SpacedSeedMaskCollection const> masks)
        : SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>(masks,
                                                                 config->genome1(),
                                                                 config->genome2(),
                                                                 idMap,
                                                                 config->matchLimit(),
                                                                 config->matchLimitDiscardSeeds(),
                                                                 config->createAllMatches(),
                                                                 config->nThreads(),
                                                                 config->quiet()),
          matches_{} {
        if (this->span() != this->weight() || this->seedSetSize() != 1) { throw std::runtime_error("[ERROR] -- SeedMapContiguous -- masks has invalid parameters"); }
    }
    void appendMatches(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> sourceMap) {
        auto sourceCast = std::static_pointer_cast<SeedMapContiguous<TwoBitKmerDataType, TwoBitSeedDataType>>(sourceMap);
        matches_.insert(matches_.end(),
                        sourceCast->matches().begin(),
                        sourceCast->matches().end());
    }
    void clearMatches() { matches_.clear(); }
    std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                   TwoBitSeedDataType>> localCopy() const;
    void createMatches() {
        this->createMatchesGeneral();
    }
    void insertMatch(KmerOccurrence const & occ1, KmerOccurrence const & occ2) {
        matches_.emplace_back(occ1, occ2);
    }
    //! Getter for member \c matches_
    auto const & matches() const { return matches_; }
    void merge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> localMap);
    //! Merge matches from a source to member \c matches_
    template <typename ContainerSource>
    void mergeMatches(ContainerSource const & matchSource) {
        this->mergeMatches_impl(matchSource, matches_);
    }
    size_t numMatches() const { return matches_.size(); }
    void quickMerge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                   TwoBitSeedDataType>> localMap);
    void output(std::ostream & outstream) const;

private:
    //! Stores k-mer matches between \c genome0 and \c genome1
    MatchesType matches_;
};

#endif // EXACTMATCHES_H
