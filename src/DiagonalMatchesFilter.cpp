#include "DiagonalMatchesFilter.h"

template <typename LinkType>
std::vector<LinkType> DiagonalMatchesFilter<LinkType>::applyDiagonalMatchesFilter(std::vector<LinkType> const & sortedMatches) {
    if (sortedMatches.size() == 0) { return std::vector<LinkType>{}; }

    // get roughly equal-sized chunks that only contain complete diagonals
    auto chunks = matchesChunks(sortedMatches);
    Timestep tsReport("Applying Diagonal Filter to Matches (" + std::to_string(chunks.size()) + " chunks)");
    std::vector<LinkType> filteredMatches;
    // spawn threads
    std::vector<std::thread> threads;
    for (auto&& chunk : chunks) {
        auto& chunkIt = chunk.first;
        auto& chunkEnd = chunk.second;
        if (chunkIt != chunkEnd) {
            threads.push_back(std::thread(&DiagonalMatchesFilter::filterMatchChunk, this,
                                          chunkIt, chunkEnd,
                                          std::ref(filteredMatches)));
        }
    }
    // wait until all threads are finished
    std::for_each(threads.begin(), threads.end(), [](std::thread & t) { t.join(); });
    tsReport.endAndPrint();
    return filteredMatches;
}



template <typename LinkType>
std::vector<std::pair<typename std::vector<LinkType>::const_iterator,
                      typename std::vector<LinkType>::const_iterator>>
  DiagonalMatchesFilter<LinkType>::matchesChunks(std::vector<LinkType> const & matches) const {
    std::vector<std::pair<typename std::vector<LinkType>::const_iterator,
                          typename std::vector<LinkType>::const_iterator>> chunks;
    auto chunkSize = static_cast<size_t>( std::ceil( static_cast<double>(matches.size()) / static_cast<double>(config_->nThreads()) ) );
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

        // advance further until chunkEnd points to new sequence pair (or past-the-end)
        while (chunkEnd != matches.end() && chunkLast->sameSequences(*chunkEnd)) {
            ++chunkEnd;
            ++chunkLast;
        }

        // chunkEnd now either a new sequence pair or past-the end
        chunks.emplace_back(matchIt, chunkEnd);
        matchIt = chunkEnd; // new chunk begin
    }

    return chunks;
}



template <typename LinkType>
void DiagonalMatchesFilter<LinkType>::filterMatchChunk(typename std::vector<LinkType>::const_iterator matchIt,
                                                       typename std::vector<LinkType>::const_iterator matchEnd,
                                                       std::vector<LinkType> & resultGlobal) {
    std::unique_lock<std::mutex> lock(mutex_, std::defer_lock);
    std::vector<LinkType> resultLocal;
    size_t skippedNotInGenome1And2Local = 0;

    while (matchIt != matchEnd) {
        auto& occ0 = matchIt->first();
        auto& occ1 = matchIt->second();
        // only interested in matches between genome0 and genome1
        if (occ0.genome() > 1 || occ1.genome() > 1) {
            ++skippedNotInGenome1And2Local;
            ++matchIt;  // move to next pair and start next iteration
            continue;
        }
        if (matchIt == matchEnd) { break; } // in case invalid match was the last one
        // if valid match, get same-sequence-pair subset
        auto firstSeqPairMatch = matchIt;
        while ((matchIt != matchEnd) && firstSeqPairMatch->sameSequences(*matchIt)) {
            ++matchIt;  // move to next pair
        }
        // process same-sequence block
        processSameSequenceMatches(firstSeqPairMatch, matchIt, resultLocal);
    }

    lock.lock();
    skippedNotInGenome1And2_ += skippedNotInGenome1And2Local;
    resultGlobal.insert(resultGlobal.end(),
                        resultLocal.begin(),
                        resultLocal.end());
    lock.unlock();
}



template <typename LinkType>
void DiagonalMatchesFilter<LinkType>::processSameSequenceMatches(typename std::vector<LinkType>::const_iterator it,
                                                                 typename std::vector<LinkType>::const_iterator end,
                                                                 std::vector<LinkType> & result) const {
    auto interSeedDistance = [](LinkType const & a, LinkType const & b) {
        auto i1= static_cast<long long>(a.first().position());
        auto i2= static_cast<long long>(b.first().position());
        auto j1= static_cast<long long>(a.second().position());
        auto j2= static_cast<long long>(b.second().position());
        auto d1 = std::abs(i2-i1);
        auto d2 = std::abs(j2-j1);
        return static_cast<size_t>(std::max(d1,d2));
    };

    if (config_->yass()) {
        // Fix previous paramters: allow overlap (true), report everything that has at least one neighbouring seed
        //   according to the new neighbouring criteria, even if they overlap by span-1
        // Start from the lowest position on the lowest diagonal, search all neighbours on close enough diagonals
        //   If neighbour proximate enough, report both
        // Should be fine not to check in both directions due to sorting, lowest diagonal groups come first,
        //   inside each group lowest positons come first so check towards lower diagonal/position neighbours
        //   was already done and reported
        std::set<LinkType> report;
        while (it != end) {
            auto neighbourCandidate = std::next(it);
            while (neighbourCandidate != end) {
                auto deltaDiagonal = neighbourCandidate->diagonal().at(1) - it->diagonal().at(1);
                auto deltaPosition = interSeedDistance(*it, *neighbourCandidate);
                if (deltaDiagonal < 0) { throw std::runtime_error("[ERROR] -- DiagonalMatchesFilter::processSameSequenceMatches -- Negative delta diagonal"); }
                if (deltaDiagonal > static_cast<long long>(config_->diagonalDelta())) {
                    break; // no use in checking the next matches as they also have too high diagonal
                } else {
                    if (deltaPosition <= config_->diagonalRho()) {
                        report.emplace(*it);
                        report.emplace(*neighbourCandidate);
                    }
                }
                ++neighbourCandidate;
            }
            ++it;
        }
        // only matches with valid neighbours in this set, perform seed grouping and report
        groupOverlappingSeeds(report);
        for (auto&& match : report) { result.emplace_back(match); }
    } else {
        // perform original M4
        std::set<LinkType> report;
        while (it != end) {
            std::vector<LinkType> sameDiag;
            sameDiag.emplace_back(*it);
            long double diagCount = 1.;
            auto neighbourCandidate = std::next(it);
            while (neighbourCandidate != end) {
                if (it->span() == 0) { throw std::runtime_error("[ERROR] -- DiagonalMatchesFilter::processSameSequenceMatches -- it span zero"); }
                if (neighbourCandidate->span() == 0) { throw std::runtime_error("[ERROR] -- DiagonalMatchesFilter::processSameSequenceMatches -- neighbourCandidate span zero"); }
                auto deltaDiagonal = neighbourCandidate->diagonal().at(1) - it->diagonal().at(1);
                if (deltaDiagonal < 0) { throw std::runtime_error("[ERROR] -- DiagonalMatchesFilter::processSameSequenceMatches -- Negative delta diagonal"); }
                if (deltaDiagonal > 0) {
                    break; // no use in checking the next matches as they also have differing diagonals
                } else {
                    auto itStart = it->position(0);
                    auto itEnd = itStart + it->span() - 1;
                    auto nStart = neighbourCandidate->position(0);
                    auto nEnd = nStart + neighbourCandidate->span() - 1;
                    bool overlap = nStart <= itEnd;
                    if (!config_->allowOverlap() && overlap) {
                        ++neighbourCandidate;
                        continue; // ignore overlapping links
                    }
                    if (nStart < itStart) { throw std::runtime_error("[ERROR] -- DiagonalMatchesFilter::processSameSequenceMatches -- Negative link distance"); }
                    auto distance = nStart - itStart;
                    if (distance <= config_->localAreaLength() && distance >= config_->minMatchDistance()) {
                        sameDiag.emplace_back(*neighbourCandidate);
                        if (overlap && (nEnd > itEnd)) {
                            // only add non-overlapping fraction to link count
                            auto len = nEnd - itEnd;
                            diagCount += static_cast<long double>(len)/static_cast<long double>(neighbourCandidate->span());
                        } else if (!overlap) {
                            diagCount += 1;
                        }
                    }
                }
                ++neighbourCandidate;
            }
            if (diagCount >= config_->diagonalThreshold()) {
                report.insert(sameDiag.begin(), sameDiag.end());
            }
            ++it;
        }
        // only matches with valid neighbours in this set, perform seed grouping and report
        groupOverlappingSeeds(report);
        for (auto&& match : report) { result.emplace_back(match); }
    }
}

template class DiagonalMatchesFilter<Link>;
template class DiagonalMatchesFilter<LinkPtr>;
