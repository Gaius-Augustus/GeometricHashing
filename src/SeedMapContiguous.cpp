#include "SeedMapContiguous.h"

using json = nlohmann::json;

template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                               TwoBitSeedDataType>> SeedMapContiguous<TwoBitKmerDataType,
                                                                      TwoBitSeedDataType>::localCopy() const {
    auto copyObj = std::make_shared<SeedMapContiguous<TwoBitKmerDataType,
                                                      TwoBitSeedDataType>>(this->span(), this->genome0_, this->genome1_,
                                                                           std::make_shared<IdentifierMapping>(*(this->idMap_)),
                                                                           this->matchLimit_, this->matchLimitDiscardSeeds_,
                                                                           this->createAllMatches_,
                                                                           this->nThreads_, this->quiet_);
    copyObj->idMap_ = this->idMap_;
    return std::static_pointer_cast<SeedMapGeneral<TwoBitKmerDataType,
                                                   TwoBitSeedDataType>>(copyObj);
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapContiguous<TwoBitKmerDataType,
                       TwoBitSeedDataType>::merge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                                                 TwoBitSeedDataType>> localMap) {
    auto localMapCast = std::static_pointer_cast<SeedMapContiguous<TwoBitKmerDataType, TwoBitSeedDataType>>(localMap);
    this->mergeGeneral(localMapCast);   // merge basic members
    matches_.insert(matches_.end(), // merge matches
                    localMapCast->matches().begin(),
                    localMapCast->matches().end());
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapContiguous<TwoBitKmerDataType,
                       TwoBitSeedDataType>::quickMerge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                                                      TwoBitSeedDataType>> localMap) {
    SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>::quickMerge(localMap);   // call base (general)
    auto localMapCast = std::static_pointer_cast<SeedMapContiguous<TwoBitKmerDataType,
                                                                   TwoBitSeedDataType>>(localMap);
    matches_.insert(matches_.end(),
                    localMapCast->matches_.begin(),
                    localMapCast->matches_.end());
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
std::pair<size_t, size_t> SeedMapContiguous<TwoBitKmerDataType,
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



template class SeedMapContiguous<TwoBitKmerDataLong, TwoBitKmerDataLong>;
template class SeedMapContiguous<TwoBitKmerDataMedium, TwoBitKmerDataMedium>;
template class SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>;
