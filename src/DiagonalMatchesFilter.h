#ifndef DIAGONALMATCHESFILTER_H
#define DIAGONALMATCHESFILTER_H

#include <algorithm>
#include <array>
#include <execution>
#include <thread>
#include <mutex>

#include "boost/dynamic_bitset.hpp"
#include "mabl3/Timestep.h"
#include <tsl/hopscotch_set.h>
#include "Configuration.h"
#include "KmerOccurrence.h"
#include "Link.h"

using namespace mabl3;

template<typename LinkType>
class DiagonalMatchesFilter {
public:
    //! c'tor
    /*! \param config Configuration object
     *
     * \details Initializes all members */
    DiagonalMatchesFilter(std::shared_ptr<Configuration const> config)
        : config_{config}, mutex_{},
          skippedNotInGenome1And2_{0} {}

    //! Apply filter to a set of matches, removing matches that didn't pass from \c matches
    std::vector<LinkType> applyDiagonalMatchesFilter(std::vector<LinkType> const & sortedLinks);
    auto config() const { return config_; }
    //! Split \c matches set in roughly equal sized chunks (one for each thread) that all contain only complete sequence pairs
    std::vector<std::pair<typename std::vector<LinkType>::const_iterator,
                          typename std::vector<LinkType>::const_iterator>> matchesChunks(std::vector<LinkType> const & matches) const;
    auto nThreads() const { return config_->nThreads(); }
    auto skipped() const { return skippedNotInGenome1And2_; }
    auto skippedNotInGenome1And2() const { return skippedNotInGenome1And2_; }
    auto span() const { return config_->span(); }

private:
    size_t processSameDiagonalMatches(std::vector<LinkType> const & sameDiagonal,
                                      boost::dynamic_bitset<> const & diagonalMatchPositions,
                                      size_t bitvectorOffset,
                                      std::vector<LinkType> & reported) const;
    void processSameSequenceMatches(typename std::vector<LinkType>::const_iterator it,
                                    typename std::vector<LinkType>::const_iterator end,
                                    std::vector<LinkType> & result) const;
    void filterMatchChunk(typename std::vector<LinkType>::const_iterator matchIt,
                          typename std::vector<LinkType>::const_iterator matchEnd,
                          std::vector<LinkType> & resultGlobal);

    std::shared_ptr<Configuration const> config_;
    std::mutex mutex_;
    size_t skippedNotInGenome1And2_;
};

#endif // DIAGONALMATCHESFILTER_H
