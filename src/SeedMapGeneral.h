#ifndef SEEDMAPGENERAL_H
#define SEEDMAPGENERAL_H

#include <cstdlib>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <type_traits>

#include "ContainerChunks.h"
#include "IdentifierMapping.h"
#include "JsonStream.h"
#include "KmerOccurrence.h"
#include "ProgressBar.h"
#include "SpacedSeedMask.h"
#include "TwoBitKmer.h"



//! Base class for a (spaced) seed map
/*! Dictates the methods a productive child needs to implement in order to work,
 * implements most methods */
template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
class SeedMapGeneral {
public:
    // outer vector: one element for each spaced seed mask
    using SeedMapType = std::unordered_map<TwoBitKmer<TwoBitSeedDataType>,
                                           std::vector<std::vector<KmerOccurrence>>,    // outer vector: one element for each mask
                                           TwoBitKmerHash<TwoBitSeedDataType>>;
    using StatisticsType = std::unordered_map<std::string,  // genome
                                              std::tuple<size_t,    // total sequence length
                                                         size_t>>;  // #non unique kmers discarded
    //! c'tor (1)
    /*! \param span k-mer length
     * \param genome0 Name of first genome
     * \param genome1 Name of second genome
     * \param globalOccurrenceThreshold Discard k-mers that globally occur more often than this
     * \param perSequenceOccurrenceThreshold Discard k-mers that occur more often than this in any sequence
     * \param createAllMatches In Triangulation Filter setting, all matches are required, so set this; otherwise can be false to save ressources
     * \param nThreads Number of threads for parallelization
     * \param quiet If true, do not print progress bars
     *
     * \details Initializes all members, creates a single SpacedSeedMask with span == weight and all 1s */
    SeedMapGeneral(size_t span, std::string const & genome0, std::string const & genome1,
                   std::shared_ptr<IdentifierMapping> idMap,
                   size_t matchLimit, bool matchLimitDiscardSeeds, bool createAllMatches,
                   size_t nThreads, bool quiet = false)
        : createAllMatches_{createAllMatches}, genome0_{genome0}, genome1_{genome1},
          idMap_{idMap}, matchLimit_{matchLimit}, matchLimitDiscardSeeds_{matchLimitDiscardSeeds},
          mutex_{}, nThreads_{nThreads}, perGenomeStats_{}, quiet_{quiet}, rd_{},
          seedMap_{}, spacedSeedMasks_{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight(span),
                                                                                  SpacedSeedMaskCollection::Span(span),
                                                                                  SpacedSeedMaskCollection::SeedSetSize(1))} {
        if (idMap_->queryGenomeID(genome0_) != 0) {
            throw std::runtime_error("[ERROR] -- SeedMapGeneral -- '--genome1' must have ID '0'");
        }
        if (idMap_->queryGenomeID(genome1_) != 1) {
            throw std::runtime_error("[ERROR] -- SeedMapGeneral -- '--genome2' must have ID '1'");
        }
    }
    //! c'tor (2)
    /*! \param masks Collection of \c seedSetSize distinct spaced seed masks of length \c span and \c weight care-positions
     * \param genome0 Name of first genome
     * \param genome1 Name of second genome
     * \param idMap Shared identifier map
     * \param globalOccurrenceThreshold Discard k-mers that globally occur more often than this
     * \param perSequenceOccurrenceThreshold Discard k-mers that occur more often than this in any sequence
     * \param createAllMatches In Triangulation Filter setting, all matches are required, so set this; otherwise can be false to save ressources
     * \param nThreads Number of threads for parallelization
     * \param quiet If true, do not print progress bars
     *
     * \details Initializes all members */
    SeedMapGeneral(std::shared_ptr<SpacedSeedMaskCollection const> masks,
                   std::string const & genome0, std::string const & genome1,
                   std::shared_ptr<IdentifierMapping> idMap,
                   size_t matchLimit, bool matchLimitDiscardSeeds, bool createAllMatches,
                   size_t nThreads, bool quiet = false)
        : createAllMatches_{createAllMatches}, genome0_{genome0}, genome1_{genome1},
          idMap_{idMap}, matchLimit_{matchLimit}, matchLimitDiscardSeeds_{matchLimitDiscardSeeds},
          mutex_{}, nThreads_{nThreads}, perGenomeStats_{}, quiet_{quiet}, rd_{},
          seedMap_{}, spacedSeedMasks_{masks} {
        if (idMap_->queryGenomeID(genome0_) != 0) {
            throw std::runtime_error("[ERROR] -- SeedMapGeneral -- '--genome1' must have ID '0'");
        }
        if (idMap_->queryGenomeID(genome1_) != 1) {
            throw std::runtime_error("[ERROR] -- SeedMapGeneral -- '--genome2' must have ID '1'");
        }
    }
    virtual ~SeedMapGeneral() {}

    //! Implementations of this method need to take care of adding a new seed to the seed map
    /*! This includes filtering and calling cleanup if neccessary, it is important that the
     * \c seedMap_ member is in a valid state according to filter criteria after calling \c addSeed() */
    virtual void addSeed(TwoBitKmer<TwoBitSeedDataType> const & seed,
                         KmerOccurrence const & occurrence,
                         size_t maskIndex);
    //! Append matches from sourceMap to this matches, children need to implement this according to their match-container
    virtual void appendMatches(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> sourceMap) = 0;
    //! Appends the json representation of a match pair to a nlohmann::json object
    void appendMatchToOutput(JsonStreamArray & jstream, KmerOccurrence const & occ0,
                             KmerOccurrence const & occ1) const;
    void applyDiagonalMatchesFilter(std::shared_ptr<std::array<size_t, 4>> skipped) {
        (void)skipped;  // intentionally unused
        throw std::runtime_error("[ERROR] -- applyDiagonalMatchesFilter does nothing in this setting, should not be called!");
    }
    //! Implementations of this method should remove all elements from the matches container member
    virtual void clearMatches() = 0;
    //! Implementations of this method should clear the seedMap member that stores the mapping from a seed to its occurrences
    virtual void clearSeedMap() { seedMap_.clear(); }
    //! Implementations of this method should create some set of seed matches across the genomes
    virtual void createMatches() = 0;
    //! Creates and inserts matches from a vector of KmerOccurrence s
    void createMatchesFromOccurrences(std::vector<KmerOccurrence> const & occurrences);
    //! Creates matches from seeds
    void createMatchesGeneral() {
        // parallel does slow things down, so use single thread
        ProgressBar pb(seedMap_.size(), quiet_);
        for (auto&& elem : seedMap_) {
            for (auto&& occurrenceVector : elem.second) {
                createMatchesFromOccurrences(occurrenceVector);
            }
            ++pb;
        }
    }
    //! Implementations of this method need to take care of adding a new seed to the seed map from a possibly longer region using a SpacedSeedMask
    /*! This includes filtering and calling cleanup if neccessary, it is important that the
     * \c seedMap_ member is in a valid state according to filter criteria after calling \c addSeed() */
    virtual void createSeed(TwoBitKmer<TwoBitKmerDataType> const & seed,
                            KmerOccurrence const & occurrence);
    //! Getter for member \c genome0_
    auto const & genome0() const { return genome0_; }
    //! Getter for member \c genome1_
    auto const & genome1() const { return genome1_; }
    //! Getter for member \c idMap_
    auto idMap() { return idMap_; }
    //! Returns a copy of \c idMap_
    auto idMapCopy() const { return *idMap_; }
    //! Insert a pairwise match, children need to implement this according to their match-container
    virtual void insertMatch(KmerOccurrence const & occ1, KmerOccurrence const & occ2) = 0;
    //! Implementations of this method should return an empty object with the same settings and blacklists but entirely distinct data
    /*! This means that e.g. shared pointer members should point to new, copied data instead of \c this data */
    virtual std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                           TwoBitSeedDataType>> localCopy() const = 0; // only exclude information copied
    //! Forward getter for member \c spacedSeedMasks_->masks()
    auto const & masks() const { return spacedSeedMasks_->masks(); }
    //! Getter for \c matchLimit_
    auto matchLimit() const { return matchLimit_; }
    //! Getter for \c matchLimitDiscardSeeds_
    auto matchLimitDiscardSeeds() const { return matchLimitDiscardSeeds_; }
    //! Implementations of this method should take care that two objects of the same type get merged into one
    /*! As with \c addSeed(), this method needs to make sure that \c seedMap_ is in a valid state afterwards */
    virtual void merge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                      TwoBitSeedDataType>> localMap) = 0;
    //! Implement merging for basic members, productive childs only need to implement merge() for their special members and call this
    /*! \c seedMap_ is in valid state after calling \c mergeGeneral() */
    void mergeGeneral(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                     TwoBitSeedDataType>> localMap);
    auto nThreads() const { return nThreads_; }
    //! Forwards to method \c numGenomes of \c idMap_
    auto numGenomes() const { return idMap_->numGenomes(); }
    //! Return number of matches
    virtual size_t numMatches() const = 0;
    //! Return number of seeds
    auto numSeeds() const { return seedMap_.size(); }
    //! Forwards to method \c numSequences of \c idMap_
    auto numSequences() const { return idMap_->numSequences(); }
    //! Implementations of this method should output all valid matches to an outstream
    virtual void output(std::ostream & outstream) const = 0;
    //! Getter for member \c perGenomeStats_
    auto const & perGenomeStatistics() const { return perGenomeStats_; }
    //! Implementations of this method should print statistics about the map creation process
    virtual void printStatistics() const;
    //! Only merge members but don't apply filtering again
    virtual void quickMerge(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                           TwoBitSeedDataType>> localMap);
    //! Getter for member \c seedSetSize_
    auto seedSetSize() const { return spacedSeedMasks_->size(); }
    //! Getter for member \c seedMap_
    auto const & seedMap() const { return seedMap_; }
    //! Getter for member spacedSeedMasks_
    auto spacedSeedMasks() const { return spacedSeedMasks_; }
    //! Getter for member \c span_
    auto span() const { return spacedSeedMasks_->maxSpan(); }
    //! Implementations of this method should update a statistics member
    virtual void updateSequenceLengthStatistics(std::string genomeName, size_t sequenceLength) {
        std::get<0>(perGenomeStats_[genomeName]) += sequenceLength;
    }
    //! Implementations of this method should update a statistics member
    virtual void updateInvalidSeedStatistics(std::string genomeName) {
        ++std::get<1>(perGenomeStats_[genomeName]);
    }
    //! Getter for member \c weight_
    auto weight() const { return spacedSeedMasks_->weight(); }

protected:
    //! Returns reference to \c seedMap_[seed], creates correctly sized vectors if new seed
    auto& accessSeedInMap(TwoBitKmer<TwoBitSeedDataType> const & seed) {
        if (seedMap_.find(seed) == seedMap_.end()) {
            seedMap_.emplace(seed, spacedSeedMasks_->size());
        }
        return seedMap_[seed];
    }
    //! Merge match containers
    template <typename ContainerSource, typename ContainerTarget>
    void mergeMatches_impl(ContainerSource const & source, ContainerTarget & target) const {
        containerChunkInsert<ContainerTarget, ContainerSource>(target, source.begin(), source.end());
    }
    //! Induces a spaced seed from a k-mer using a SpacedSeedMask
    TwoBitKmer<TwoBitSeedDataType> spacedSeedFromKmer(TwoBitKmer<TwoBitKmerDataType> const & kmer,
                                                      SpacedSeedMask const & mask) const;
    //! Induces a spaced seed from a k-mer using a SpacedSeedMask via its index
    TwoBitKmer<TwoBitSeedDataType> spacedSeedFromKmer(TwoBitKmer<TwoBitKmerDataType> const & kmer,
                                                      size_t maskIndex) const {
        return spacedSeedFromKmer(kmer, spacedSeedMasks_->masks().at(maskIndex));
    }

    //! Flag if all matches should be created, if false, only create matches from seeds that occur in both genome 0 and 1
    bool createAllMatches_;
    //! First genome from which matches are sought
    std::string const genome0_;
    //! Second genome from which matches are sought
    std::string const genome1_;
    //! Identifier Mapping
    std::shared_ptr<IdentifierMapping> idMap_;
    //! If number of possible matches from a seed exceeds this number, sample this many of all possible matches
    size_t matchLimit_;
    //! Discard seeds exceeding \c matchLimit_ rather than sampling
    bool matchLimitDiscardSeeds_;
    //! General lock
    std::mutex mutex_;
    //! Number of threads for parallel processing
    size_t const nThreads_;
    //! Stores statistics about the seedMap creation
    StatisticsType perGenomeStats_;
    //! Don't show progress bars
    bool quiet_;
    //! Used to obtain seed for random link selection if there are too many possibilities
    std::random_device rd_;
    //! seedMap
    SeedMapType seedMap_;
    //! Vector of all SpacedSeedMask s
    std::shared_ptr<SpacedSeedMaskCollection const> spacedSeedMasks_;
};

#endif // SEEDMAPGENERAL_H
