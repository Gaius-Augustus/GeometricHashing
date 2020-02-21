#ifndef DIAGONALMATCHESFILTER_H
#define DIAGONALMATCHESFILTER_H

#include <array>
#include <thread>
#include <mutex>

#include "boost/dynamic_bitset.hpp"
#include "Configuration.h"
#include "KmerOccurrence.h"
#include "Link.h"
#include "Timestep.h"

template<typename KmerOccurrencePairType>
class DiagonalMatchesFilter {
public:
    //! c'tor (1)
    /*! \param config Configuration object
     *
     * \details Initializes all members */
    DiagonalMatchesFilter(std::shared_ptr<Configuration const> config)
        : allowOverlap_{config->allowOverlap()},
          diagonalThreshold_{config->diagonalThreshold()},
          localAreaLength_{config->localAreaLength()},
          minMatchDistance_{config->minMatchDistance()},
          mutex_{}, nThreads_{config->nThreads()}, quiet_{config->quiet()},
          skippedNotInGenome1And2_{0}, skippedOverlappedOrTooClose_{0},
          skippedTooFewDiagonalElements_{0}, skippedTooFewNeighbours_{0},
          span_{config->span()} {}
    //! c'tor (2)
    /*! \param localAreaLength Search in this area around a match for neighbouring matches
     * \param minMatchDistance When looking for neighbouring matches, discard all matches that are closer than this to the match
     * \param allowOverlap Disables \c minMatchDistance, allow neighbouring matches to overlap with the match, counting only no-overlapping fractions
     * \param nThreads Apply filter in parallel on this many threads
     * \param quiet Do not show progress bars if \c true
     *
     * \details Initializes all members */
    DiagonalMatchesFilter(size_t span, double diagonalThreshold,
                          size_t localAreaLength, size_t minMatchDistance,
                          bool allowOverlap, size_t nThreads, bool quiet = false)
        : allowOverlap_{allowOverlap}, diagonalThreshold_{diagonalThreshold},
          localAreaLength_{localAreaLength}, minMatchDistance_{minMatchDistance},
          mutex_{}, nThreads_{nThreads}, quiet_{quiet},
          skippedNotInGenome1And2_{0}, skippedOverlappedOrTooClose_{0},
          skippedTooFewDiagonalElements_{0}, skippedTooFewNeighbours_{0},
          span_{span} {}

    auto allowOverlap() const { return allowOverlap_; }
    //! Apply filter to an already ordered set of matches
    /*! It must be possible to get \c const_iterators from \c matches
     * It also must be possible to get a KmerOccurrenceDistance object from each element in \c matches */
    void applyDiagonalMatchesFilter(std::set<KmerOccurrencePairType> const & matches,
                                    std::vector<KmerOccurrencePairType> & result,
                                    bool quiet = false);
    auto diagonalThreshold() const{ return diagonalThreshold_; }
    auto localAreaLength() const { return localAreaLength_; }
    //! Split \c matches set in roughly equal sized chunks (one for each thread) that all contain only complete diagonals
    std::vector<std::pair<typename std::set<KmerOccurrencePairType>::const_iterator,
                          typename std::set<KmerOccurrencePairType>::const_iterator>> matchesChunks(std::set<KmerOccurrencePairType> const & matches) const;
    auto minMatchDistance() const { return minMatchDistance_; }
    auto nThreads() const { return nThreads_; }
    auto skipped() const { return (skippedNotInGenome1And2_ + skippedOverlappedOrTooClose_
                                   + skippedTooFewDiagonalElements_ + skippedTooFewNeighbours_); }
    auto skippedNotInGenome1And2() const { return skippedNotInGenome1And2_; }
    auto skippedOverlappedOrTooClose() const { return skippedOverlappedOrTooClose_; }
    auto skippedTooFewDiagonalElements() const { return skippedTooFewDiagonalElements_; }
    auto skippedTooFewNeighbours() const { return skippedTooFewNeighbours_; }
    auto span() const { return span_; }

private:
    size_t processSameDiagonalMatches(std::vector<KmerOccurrencePairType> const & sameDiagonal,
                                      boost::dynamic_bitset<> const & diagonalMatchPositions,
                                      size_t bitvectorOffset,
                                      std::vector<KmerOccurrencePairType> & reported) const;
    void reportMatchChunk(typename std::set<KmerOccurrencePairType>::const_iterator matchIt,
                          typename std::set<KmerOccurrencePairType>::const_iterator matchEnd,
                          std::vector<KmerOccurrencePairType> & reportedGlobal);

    bool allowOverlap_;
    double diagonalThreshold_;
    size_t localAreaLength_;
    size_t minMatchDistance_;
    std::mutex mutex_;
    size_t nThreads_;
    bool quiet_;
    size_t skippedNotInGenome1And2_;
    size_t skippedOverlappedOrTooClose_;
    size_t skippedTooFewDiagonalElements_;
    size_t skippedTooFewNeighbours_;
    size_t span_;
};

#endif // DIAGONALMATCHESFILTER_H
