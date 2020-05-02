#include "SeedMapSpaced.h"

using json = nlohmann::json;

template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapSpaced<TwoBitKmerDataType,
                   TwoBitSeedDataType>::applyDiagonalMatchesFilter(std::shared_ptr<std::array<size_t, 4>> skipped) {
    std::vector<KmerOccurrencePair> reported;
    auto filter = DiagonalMatchesFilter<KmerOccurrencePair>(config_);
    std::vector<KmerOccurrencePair> orderedMatches;
    Timestep tsOrdered("Creating ordered list of matches");
    orderedMatches.insert(orderedMatches.begin(), matches_.begin(), matches_.end());
    matches_.clear();
    Timestep tsSort("Sorting");
    std::sort(std::execution::par, orderedMatches.begin(), orderedMatches.end());
    tsSort.endAndPrint();
    tsOrdered.endAndPrint();

    Timestep ts("Filtering");
    filter.applyDiagonalMatchesFilter(orderedMatches, reported);
    (*skipped)[0] = filter.skippedNotInGenome1And2();
    (*skipped)[1] = filter.skippedOverlappedOrTooClose();
    (*skipped)[2] = filter.skippedTooFewDiagonalElements();
    (*skipped)[3] = filter.skippedTooFewNeighbours();

    std::cout << "[INFO] -- SeedMapSpaced::reportMatches -- Skipped " << filter.skippedNotInGenome1And2() << " matches that not included genome1 and 2";
    std::cout << ", " << filter.skippedOverlappedOrTooClose() << " matches that overlapped or were too close to another match";
    std::cout << ", " << filter.skippedTooFewDiagonalElements() << " matches from diagonals with less than --diagonal-threshold = " << filter.diagonalThreshold() << " matches and ";
    std::cout << filter.skippedTooFewNeighbours() << " matches that did not have enough neighbouring matches ";
    std::cout << "(total: skipped " << filter.skipped() << " / " << matches_.size() << " matches)" << std::endl << std::endl;

    orderedMatches.clear();
    matches_.insert(reported.begin(), reported.end());
    ts.endAndPrint();
}



template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                               TwoBitSeedDataType>> SeedMapSpaced<TwoBitKmerDataType,
                                                                  TwoBitSeedDataType>::localCopy() const {
    auto copyObj = std::make_shared<SeedMapSpaced<TwoBitKmerDataType,
                                                  TwoBitSeedDataType>>(std::make_shared<Configuration const>(*config_),
                                                                       std::make_shared<IdentifierMapping>(*(this->idMap_)),
                                                                       std::make_shared<SpacedSeedMaskCollection const>(*(this->spacedSeedMasks_)));
    return std::static_pointer_cast<SeedMapGeneral<TwoBitKmerDataType,
                                                   TwoBitSeedDataType>>(copyObj);
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapSpaced<TwoBitKmerDataType,
                   TwoBitSeedDataType>::merge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                                             TwoBitSeedDataType>> localMap) {
    auto localMapCast = std::static_pointer_cast<SeedMapSpaced<TwoBitKmerDataType, TwoBitSeedDataType>>(localMap);
    this->mergeGeneral(localMapCast);   // merge basic members
    // merge matches
    matches_.insert(localMapCast->matches().begin(),
                    localMapCast->matches().end());
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapSpaced<TwoBitKmerDataType,
                   TwoBitSeedDataType>::quickMerge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                                                  TwoBitSeedDataType>> localMap) {
    SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>::quickMerge(localMap);   // call base (general)
    auto localMapCast = std::static_pointer_cast<SeedMapSpaced<TwoBitKmerDataType,
                                                               TwoBitSeedDataType>>(localMap);
    matches_.insert(localMapCast->matches_.begin(),
                    localMapCast->matches_.end());
}



template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
std::pair<size_t, size_t> SeedMapSpaced<TwoBitKmerDataType,
                                        TwoBitSeedDataType>::output(std::ostream & outstream) const {
    auto jstream = JsonStreamArray(outstream);
    size_t skipped = 0;
    size_t written = 0;
    for (auto&& elem : matches_) {
        auto& occ0 = elem.first();
        auto& occ1 = elem.second();
        // only interested in matches between genome0 and genome1
        if (occ0.genome() > 1 || occ1.genome() > 1) {
            ++skipped;
            continue;
        }
        this->appendMatchToOutput(jstream, occ0, occ1);
        ++written;
    }
    return std::pair<size_t, size_t>{written, skipped};
}



template class SeedMapSpaced<TwoBitKmerDataLong, TwoBitKmerDataLong>;
template class SeedMapSpaced<TwoBitKmerDataLong, TwoBitKmerDataMedium>;
template class SeedMapSpaced<TwoBitKmerDataLong, TwoBitKmerDataShort>;
template class SeedMapSpaced<TwoBitKmerDataMedium, TwoBitKmerDataMedium>;
template class SeedMapSpaced<TwoBitKmerDataMedium, TwoBitKmerDataShort>;
template class SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>;
