#include "DiagonalMatchesFilter.h"

template<typename KmerOccurrencePairType>
void DiagonalMatchesFilter<KmerOccurrencePairType>::applyDiagonalMatchesFilter(std::vector<KmerOccurrencePairType> const & matches,
                                                                               std::vector<KmerOccurrencePairType> & result,
                                                                               bool quiet) {
    if (matches.size() == 0) { return; }

    // get roughly equal-sized chunks that only contain complete diagonals
    auto chunks = matchesChunks(matches);

    std::unique_ptr<Timestep> tsReport;
    if (!quiet) {
      tsReport = std::make_unique<Timestep>("Applying Diagonal Filter to Matches (" + std::to_string(chunks.size()) + " chunks)");
    }

    // spawn threads
    std::vector<std::thread> threads;
    for (auto&& chunk : chunks) {
        auto& chunkIt = chunk.first;
        auto& chunkEnd = chunk.second;
        if (chunkIt != chunkEnd) {
            threads.push_back(std::thread(&DiagonalMatchesFilter<KmerOccurrencePairType>::reportMatchChunk, this,
                                          chunkIt, chunkEnd,
                                          std::ref(result)));
        }
    }

    // wait until all threads are finished
    std::for_each(threads.begin(), threads.end(), [](std::thread & t) { t.join(); });
    if (!quiet) {
        tsReport->endAndPrint(Timestep::seconds);
    }
    return;
}



template<typename KmerOccurrencePairType>
std::vector<std::pair<typename std::vector<KmerOccurrencePairType>::const_iterator,
                      typename std::vector<KmerOccurrencePairType>::const_iterator>>
  DiagonalMatchesFilter<KmerOccurrencePairType>::matchesChunks(std::vector<KmerOccurrencePairType> const & matches) const {
    std::vector<std::pair<typename std::vector<KmerOccurrencePairType>::const_iterator,
                          typename std::vector<KmerOccurrencePairType>::const_iterator>> chunks;
    auto chunkSize = static_cast<size_t>( std::ceil( static_cast<double>(matches.size()) / static_cast<double>(nThreads_) ) );
    auto matchIt = matches.begin();
    auto chunkLast = matches.begin();
    auto chunkEnd = matches.begin();
    ++chunkEnd;

    // special case: only one match
    if (matches.size() < 2) {
        chunks.emplace_back(matches.begin(), matches.end());
        return chunks;
    }

    while (chunkEnd != matches.end()) { // if matches.size() == 1, this is not executed!
        // advance new chunk last and end until chunk size or matches end
        for (size_t i = 0; i < chunkSize; ++i) {
            if (chunkEnd == matches.end()) { break; }
            ++chunkEnd;
            ++chunkLast;
        }

        // advance further until chunkEnd points to new diagonal (or past-the-end)
        while (chunkEnd != matches.end() && chunkLast->sameDistance(*chunkEnd)) {
            ++chunkEnd;
            ++chunkLast;
        }

        // chunkEnd now either a new diagonal or past-the end
        chunks.emplace_back(matchIt, chunkEnd);
        matchIt = chunkEnd; // new chunk begin
    }

    return chunks;
}



template<typename KmerOccurrencePairType>
void DiagonalMatchesFilter<KmerOccurrencePairType>::reportMatchChunk(typename std::vector<KmerOccurrencePairType>::const_iterator matchIt,
                                                                     typename std::vector<KmerOccurrencePairType>::const_iterator matchEnd,
                                                                     std::vector<KmerOccurrencePairType> & reportedGlobal) {
    std::unique_lock<std::mutex> lock(mutex_, std::defer_lock);
    std::vector<KmerOccurrencePairType> reportedLocal;
    size_t skippedNotInGenome1And2Local = 0;
    size_t skippedOverlappedOrTooCloseLocal = 0;
    size_t skippedTooFewDiagonalElementsLocal = 0;
    size_t skippedTooFewNeighboursLocal = 0;
    std::vector<KmerOccurrencePairType> sameDiagonal;

    while (matchIt != matchEnd) {
        auto& occ0 = matchIt->first();
        auto& occ1 = matchIt->second();
        // only interested in matches between genome0 and genome1
        if (occ0.genome() > 1 || occ1.genome() > 1) {
            ++skippedNotInGenome1And2Local;
            ++matchIt;  // move to next pair and start next iteration
            continue;
        }
        // if valid match, get same-diagonal subsets (are automatically all valid)
        auto firstDiagonalMatch = *matchIt;
        size_t nextAllowedPosition = 0;
        while ((matchIt != matchEnd) && firstDiagonalMatch.sameDistance(*matchIt)) {
            if (allowOverlap_ || matchIt->first().position() >= nextAllowedPosition) {  // filter out overlapping matches if appropriate
                sameDiagonal.emplace_back(*matchIt);
                nextAllowedPosition = matchIt->first().position() + this->span_ + minMatchDistance_;
            } else {
                ++skippedOverlappedOrTooCloseLocal;
            }
            ++matchIt;  // move to next pair
        }
        // sameDiagonal now contains all elements on this respective diagonal
        // matchIt points to the first element of a different diagonal (or to past-the-end)
        if (static_cast<double>(sameDiagonal.size()) < diagonalThreshold_) {
            skippedTooFewDiagonalElementsLocal += sameDiagonal.size();
            sameDiagonal.clear();
            continue;   // not enough matches, next diagonal
        } else {
            auto diagonalLength = sameDiagonal.rbegin()->first().position() + this->span_ - sameDiagonal.begin()->first().position();
            auto bitvectorOffset = sameDiagonal.begin()->first().position();
            boost::dynamic_bitset<> diagonalMatchPositions(diagonalLength);
            for (auto&& match : sameDiagonal) {
                diagonalMatchPositions.set(match.first().position() - bitvectorOffset, this->span_, true);
            }
            // check for each match if enough valid neighbouring matches exist
            skippedTooFewNeighboursLocal += processSameDiagonalMatches(sameDiagonal,
                                                                       diagonalMatchPositions,
                                                                       bitvectorOffset,
                                                                       reportedLocal);
            sameDiagonal.clear();
        }
        // matchIt already points to next diagonal (or past-the-end)
    }

    lock.lock();
    skippedNotInGenome1And2_ += skippedNotInGenome1And2Local;
    skippedOverlappedOrTooClose_ += skippedOverlappedOrTooCloseLocal;
    skippedTooFewDiagonalElements_ += skippedTooFewDiagonalElementsLocal;
    skippedTooFewNeighbours_ += skippedTooFewNeighboursLocal;
    reportedGlobal.insert(reportedGlobal.end(),
                          reportedLocal.begin(),
                          reportedLocal.end());
    lock.unlock();
}



template<typename KmerOccurrencePairType>
size_t DiagonalMatchesFilter<KmerOccurrencePairType>::processSameDiagonalMatches(std::vector<KmerOccurrencePairType> const & sameDiagonal,
                                                                                 boost::dynamic_bitset<> const & diagonalMatchPositions,
                                                                                 size_t bitvectorOffset,
                                                                                 std::vector<KmerOccurrencePairType> & reported) const {
    size_t skipped = 0;
    auto halfArea = static_cast<size_t>(std::ceil(static_cast<double>(localAreaLength_)/2.));
    for (auto&& match : sameDiagonal) {
        auto center = KmerOccurrence::centerPosition(match.first().position(), this->span_) - bitvectorOffset;

        auto leftBorder = (halfArea >= center) ? 0 : center - halfArea;
        auto rightBorder = ((halfArea + center) >= diagonalMatchPositions.size())
                            ? diagonalMatchPositions.size() - 1
                            : center + halfArea;
        size_t bitcount = 0;
        for (auto i = leftBorder; i <= rightBorder; ++i) {
            if (diagonalMatchPositions.test(i)) { ++bitcount; }
        }
        size_t seedCount = std::floor(static_cast<double>(bitcount)/static_cast<double>(this->span_));
        if (seedCount >= diagonalThreshold_) {
            reported.emplace_back(match);
        } else {
            ++skipped;
        }
    }
    return skipped;
}



template class DiagonalMatchesFilter<KmerOccurrencePair>;
template class DiagonalMatchesFilter<Link>;
