#include "SeedMapGeneral.h"

template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
TwoBitKmer<TwoBitSeedDataType> SeedMapGeneral<TwoBitKmerDataType,
                                              TwoBitSeedDataType>::spacedSeedFromKmer(TwoBitKmer<TwoBitKmerDataType> const & kmer,
                                                                                      SpacedSeedMask const & mask) const {
    auto setBits = mask.getSetPositions();
    std::string inducedKmer;
    for (auto&& pos : setBits) {
        inducedKmer += kmer.base(pos);
    }
    return TwoBitKmer<TwoBitSeedDataType>(inducedKmer);
}



template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapGeneral<TwoBitKmerDataType,
                    TwoBitSeedDataType>::addSeed(TwoBitKmer<TwoBitSeedDataType> const & seed,
                                                 KmerOccurrence const & occurrence,
                                                 size_t maskIndex) {
    if (seed.length() != spacedSeedMasks_->weight()) { throw std::runtime_error("[ERROR] -- SeedMapGeneral::addSeed() -- Wrong seed length"); }
    this->accessSeedInMap(seed).at(maskIndex).emplace_back(occurrence);
}



template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapGeneral<TwoBitKmerDataType,
                    TwoBitSeedDataType>::createSeed(TwoBitKmer<TwoBitKmerDataType> const & seed,
                                                    KmerOccurrence const & occurrence) {
    if (seed.length() != this->span()) { throw std::runtime_error("[ERROR] -- SeedMapGeneral::createSeed() -- Wrong seed length"); }
    for (size_t i = 0; i < spacedSeedMasks_->size(); ++i) {
        addSeed(spacedSeedFromKmer(seed, spacedSeedMasks_->masks()[i]), occurrence, i);
    }
}



template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapGeneral<TwoBitKmerDataType,
                    TwoBitSeedDataType>::createMatchesFromOccurrences(std::vector<KmerOccurrence> const & occurrences) {
    auto calculatePossibleMatches = [](std::vector<KmerOccurrence> const & occurrences) {
        std::map<size_t, size_t> occurrencesPerGenome;
        // check how often occurred in which genome
        for (auto&& occurrence : occurrences) {
            if (occurrencesPerGenome.find(occurrence.genome()) == occurrencesPerGenome.end()) {
                occurrencesPerGenome[occurrence.genome()] = 1;
            } else {
                ++occurrencesPerGenome[occurrence.genome()];
            }
        }
        size_t nmatches = 0;
        for (auto it = occurrencesPerGenome.begin(); it != occurrencesPerGenome.end(); ++it) {
            auto it2 =it;
            ++it2;
            for (; it2 != occurrencesPerGenome.end(); ++it2) {
                nmatches += (*it).second * (*it2).second;
            }
        }
        return nmatches;
    };

    auto matchmaker = [this](std::vector<KmerOccurrence> const & occs, size_t nPossibleMatches) {
        // classically, create all possible matches
        if (matchLimit_ >= nPossibleMatches) {
            for (auto it1 = occs.begin(); it1 != occs.end(); ++it1) {
                auto it2 = it1;
                for (++it2; it2 != occs.end(); ++it2) {
                    if (it1->genome() != it2->genome()) {
                        insertMatch(*it1, *it2);
                    }
                }
            }
        } else if (!matchLimitDiscardSeeds_) {
            // otherwise, sample until enough matches were created
            std::set<std::pair<KmerOccurrence, KmerOccurrence>> seenMatches;
            std::vector<KmerOccurrence> pairCandidate;
            size_t insertedMatches = 0;
            std::default_random_engine rng(rd_()); // create default rng with seed from rd_
            while (insertedMatches < matchLimit_) {
                std::sample(occs.begin(), occs.end(),
                            std::back_inserter(pairCandidate), 2, rng);
                if (pairCandidate.at(0).genome() != pairCandidate.at(1).genome()) {
                    std::sort(pairCandidate.begin(), pairCandidate.end());
                    auto pair = std::pair<KmerOccurrence, KmerOccurrence>{pairCandidate.at(0), pairCandidate.at(1)};
                    if (seenMatches.find(pair) == seenMatches.end()) {
                        insertMatch(pair.first, pair.second);
                        seenMatches.emplace(pair);
                        ++insertedMatches;
                    }
                }
                pairCandidate.clear();
            }
        }
    };

    if (occurrences.size() >= 2) { // do nothing if only one occurrence
        if (!createAllMatches_ && idMap_->numGenomes() > 2) {
            std::vector<KmerOccurrence> relevantOccurrences;
            for (auto&& occurrence : occurrences) {
                if (occurrence.genome() <= 1) { relevantOccurrences.emplace_back(occurrence); }
            }
            auto nPossibleMatches = calculatePossibleMatches(relevantOccurrences);
            matchmaker(relevantOccurrences, nPossibleMatches);
        } else {
            auto nPossibleMatches = calculatePossibleMatches(occurrences);
            matchmaker(occurrences, nPossibleMatches);
        }
    }
}



template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapGeneral<TwoBitKmerDataType,
                    TwoBitSeedDataType>::appendMatchToOutput(JsonStreamArray & jstream,
                                                             KmerOccurrence const & occ0,
                                                             KmerOccurrence const & occ1) const {
    std::array<std::string, 4> occArray0{
        std::to_string(KmerOccurrence::centerPosition(occ0.position(), this->span())),
        std::to_string(occ0.reverse()),
        this->idMap_->queryGenomeName(occ0.genome()),
        this->idMap_->querySequenceName(occ0.sequence())
    };
    std::array<std::string, 4> occArray1 {
        std::to_string(KmerOccurrence::centerPosition(occ1.position(), this->span())),
        std::to_string(occ1.reverse()),
        this->idMap_->queryGenomeName(occ1.genome()),
        this->idMap_->querySequenceName(occ1.sequence())
    };
    jstream << std::array<std::array<std::string, 4>, 2>{occArray0, occArray1};
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapGeneral<TwoBitKmerDataType,
                    TwoBitSeedDataType>::mergeGeneral(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                                                     TwoBitSeedDataType>> localMap) {
    auto localMapCast = std::static_pointer_cast<SeedMapGeneral<TwoBitKmerDataType,
                                                                TwoBitSeedDataType>>(localMap);
    // merge perGenomeStats_
    for (auto&& node : localMapCast->perGenomeStats_) {
        if (this->perGenomeStats_.find(node.first) != this->perGenomeStats_.end()) {
            std::get<0>(this->perGenomeStats_[node.first]) += std::get<0>(node.second);
            std::get<1>(this->perGenomeStats_[node.first]) += std::get<1>(node.second);
        } else {
            this->perGenomeStats_[node.first] = node.second;
        }
    }
    // merge seedMap_
    for (auto&& node : localMapCast->seedMap_) {
        auto& seed = node.first;
        for (size_t i = 0; i < spacedSeedMasks_->size(); ++i) {
            auto& localOccurrenceVector = node.second.at(i);
            if (this->seedMap_.find(seed) == this->seedMap_.end()) {
                this->accessSeedInMap(seed).at(i) = localOccurrenceVector;
            } else {
                auto& globalOccurrenceVector = this->seedMap_[seed].at(i);
                // each local instance starts empty and gets entirely different sequences,
                // so no overlaps and thus no need for set_union
                globalOccurrenceVector.insert(globalOccurrenceVector.end(),
                                              localOccurrenceVector.begin(), localOccurrenceVector.end());
            }
        }
    }
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapGeneral<TwoBitKmerDataType,
                    TwoBitSeedDataType>::quickMerge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                                                   TwoBitSeedDataType>> localMap) {
    // merge perGenomeStats_
    for (auto&& node : localMap->perGenomeStats_) {
        if (this->perGenomeStats_.find(node.first) != this->perGenomeStats_.end()) {
            std::get<0>(this->perGenomeStats_[node.first]) += std::get<0>(node.second);
            std::get<1>(this->perGenomeStats_[node.first]) += std::get<1>(node.second);
        } else {
            this->perGenomeStats_[node.first] = node.second;
        }
    }
    // merge seedMap_
    for (auto&& elem : localMap->seedMap_) {
        auto& seed = elem.first;
        for (size_t i = 0; i < spacedSeedMasks_->size(); ++i) {
            auto& localOccurrenceVector = elem.second.at(i);
            this->accessSeedInMap(seed).at(i).insert(this->seedMap_.at(seed).at(i).end(),
                                                     localOccurrenceVector.begin(), localOccurrenceVector.end());
        }
    }
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void SeedMapGeneral<TwoBitKmerDataType,
                    TwoBitSeedDataType>::printStatistics() const {
    auto numKmers = seedMap_.size();
    std::cout << "Number of unique k-mers created from the input files: " << numKmers << std::endl;
    for (auto&& genomeStat : perGenomeStats_) {
        auto totalSequenceLength = std::get<0>(genomeStat.second);
        auto nonUniqueKmersDiscarded = std::get<1>(genomeStat.second);
        auto approxPercent = (totalSequenceLength == 0)
                ? 0
                : std::ceil(100.*(static_cast<double>(nonUniqueKmersDiscarded)/static_cast<double>(totalSequenceLength)));
        std::cout << "Genome " << genomeStat.first << ":" << std::endl;
        std::cout << "\tTotal sequence length (sum of all sequence lengths): " << totalSequenceLength << std::endl;
        std::cout << "\t" << nonUniqueKmersDiscarded << " non-unique k-mers not read from the sequences, i.e." << std::endl;
        std::cout << "\tapprox. " << approxPercent << "% of the sequences discarded" << std::endl;
    }
}


template class SeedMapGeneral<TwoBitKmerDataLong, TwoBitKmerDataLong>;
template class SeedMapGeneral<TwoBitKmerDataLong, TwoBitKmerDataMedium>;
template class SeedMapGeneral<TwoBitKmerDataLong, TwoBitKmerDataShort>;
template class SeedMapGeneral<TwoBitKmerDataMedium, TwoBitKmerDataMedium>;
template class SeedMapGeneral<TwoBitKmerDataMedium, TwoBitKmerDataShort>;
template class SeedMapGeneral<TwoBitKmerDataShort, TwoBitKmerDataShort>;
