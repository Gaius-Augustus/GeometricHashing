#ifndef EXTRACTSEEDS_H
#define EXTRACTSEEDS_H

#include <filesystem>
#include <functional>
#include <iterator>
#include <memory>
#include <mutex>
#include <thread>

#include "mabl3/ProgressBar.h"
#include "mabl3/Timestep.h"
#include "prettyprint.hpp"
#include "Configuration.h"
#include "ContainerChunks.h"
#include "FastaCollection.h"
#include "FastaRepresentation.h"
#include "IdentifierMapping.h"
#include "KmerOccurrence.h"
#include "ParallelizationUtils.h"
#include "ParallelProgressBarHandler.h"
#include "TwoBitKmer.h"

using namespace mabl3;



//! Reads sequences and extracts seeds from it, feeding seeds into SeedMap
template <typename TwoBitSeedDataType>
class ExtractSeeds {
public:
    //! c'tor (1)
    /*! \param fastaCollection FastaCollection of input fasta files
     * \param seedMap Pointer to SeedMap
     * \param nThreads Number of threads to use for parallel sequence processing
     * \param noProgressbar Do not print progress bars
     * \param skipMatchCreation Do not create matches after extracting seeds from sequences
     *
     * \details Reads one fasta file at a time, extracting seeds from its sequences in parallel and feeding it into \c seedMap.
     * Calls \c createMatches() from \c seedMap so results are ready after this is done */
    ExtractSeeds(std::shared_ptr<SpacedSeedMaskCollection const> maskCollection,
                 std::shared_ptr<IdentifierMapping const> idMap,
                 std::shared_ptr<Configuration const> config)
        : config_{config}, idMap_{idMap},
          maskCollection_{maskCollection}, mutexMemberAccess_{},
          mutexOutput_{}, stringhasher_{} {}

    //! Create all spaced seeds from a valid k-mer and call the callback for inserting that seed into a SeedMapType
    template <typename SeedInsertCallbackType>
    void createSeed(std::string const & kmer, SeedInsertCallbackType insertValue,
                    std::function<void(TwoBitKmer<TwoBitSeedDataType>, SeedInsertCallbackType, size_t)> const & seedInsertCallback) const {
        if (kmer.length() < maskCollection_->maxSpan()) { throw std::runtime_error("[ERROR] -- ExtractSeeds::createSeed() -- Kmer too short"); }
        for (char c : kmer) {
            if ((c != 'A') && (c != 'C') && (c != 'G') && (c != 'T')) { return; }
        }
        for (size_t i = 0; i < maskCollection_->size(); ++i) {
            // create spaced seed
            bool lowComplexity = false;
            auto seed = seedFromKmer(kmer, maskCollection_->masks().at(i), lowComplexity);
            // insert seed, possibly checking for low complexity
            if ((!lowComplexity) || (!config_->redmask())) { seedInsertCallback(seed, insertValue, i); }
        }
    }
    //! Run seed extraction from fasta
    void extractFromFastas(std::shared_ptr<FastaCollection const> fastaCollection,
                           std::function<void(TwoBitKmer<TwoBitSeedDataType>, KmerOccurrence, size_t)> seedInsertCallback,
                           bool parallel = true) {
        bool silent = (!parallel) || (config_->verbose() == 0); // assume that this is called from multiple threads if parallel not allowed, thus run silently
        // Search for genome0 and 1 in inputFiles and abort if not found
        validateInputFiles(*fastaCollection);
        // Run extraction
        auto wrapper = [this,
                        &seedInsertCallback](ParallelProgressBar & pb,
                                             typename FastaRepresentation::FastaMapType::const_iterator sequencesIt,
                                             typename FastaRepresentation::FastaMapType::const_iterator sequencesEnd){
            readSequenceParallel(seedInsertCallback, pb, sequencesIt, sequencesEnd);
        };
        // read files one by one and fill seedMap_
        for (auto&& elem : fastaCollection->collection()) {
            auto& fasta = elem.second;
            size_t nthreads = parallel ? config_->nThreads() : 1;
            Timestep ts("Extracting k-mers from " + fasta.genome() + " on "
                        + std::to_string(nthreads) + " threads",
                        silent);
            ParallelProgressBar pb(fasta.numSequences(), (silent || (config_->verbose() < 2)));
            if (parallel && config_->nThreads() > 1) {
                executeParallel(fasta.headerToSequence(),
                                config_->nThreads(),
                                wrapper, std::ref(pb));
            } else {
                readSequences(fasta, seedInsertCallback, pb);
            }
            pb.unprotectedProgressBar().finish();
            ts.endAndPrint();
        }
    }
    //! Run seed extraction from metagraph
    void extractFromMetagraph(std::function<void(TwoBitKmer<TwoBitSeedDataType>, MetagraphInterface::NodeID, size_t)> seedInsertCallback) {
        ParallelProgressBar pb{config_->metagraphInterface()->numNodes(), config_->verbose() < 2};
        auto callback = [this, &seedInsertCallback, &pb](std::string const & kmer, MetagraphInterface::NodeID nodeID) {
            if (discardKmer(kmer)) { return; }
            createSeed<MetagraphInterface::NodeID>(kmer, nodeID, seedInsertCallback);
            if (!(config_->verbose() == 2)) { pb.increase(); }
        };
        config_->metagraphInterface()->iterateNodes(callback);
        pb.unprotectedProgressBar().finish();
    }
    static TwoBitKmer<TwoBitSeedDataType> seedFromKmer(std::string const & kmer, SpacedSeedMask const & mask,
                                                       bool & lowComplexity) {
        tsl::hopscotch_set<char> bases;
        auto setBits = mask.getSetPositions();
        std::string inducedSeed;
        for (auto&& pos : setBits) {
            inducedSeed += kmer.at(pos);
            bases.emplace(kmer.at(pos));
        }
        lowComplexity = (bases.size() <= 2);
        return TwoBitKmer<TwoBitSeedDataType>{inducedSeed};
    }

protected:
    //! Check if thinning leads to discarding this kmer
    bool discardKmer(std::string const & kmer) const {
        auto hash = stringhasher_(kmer);
        return ((hash % config_->thinning()) == 1);
    }
    //! Extract seeds from a single sequence
    void extractSeeds(std::string const & sequence,
                      size_t sequenceID,
                      size_t genomeID,
                      std::string const & sequenceName,
                      std::function<void(TwoBitKmer<TwoBitSeedDataType>, KmerOccurrence, size_t)> seedInsertCallback,
                      std::unique_lock<std::mutex> & outputLock) {
        if (sequence.size() < maskCollection_->maxSpan()) {
            outputLock.lock();
            std::cerr << "[WARNING] -- ExactKmerMap -- Sequence " << sequenceName
                      << " is shorter than l (" << maskCollection_->maxSpan() << ")" << std::endl;
            outputLock.unlock();
        } else {
            for (size_t p = 0; p < (sequence.size() - maskCollection_->maxSpan() + 1); ++p) {
                auto kmer = sequence.substr(p, maskCollection_->maxSpan());
                if (discardKmer(kmer)) { continue; }
                //TwoBitKmer<TwoBitKmerDataType> twoBitKmer(kmer, validFlag);
                auto occurrence = KmerOccurrence(genomeID, sequenceID, p, false, kmer);
                createSeed(kmer, occurrence, seedInsertCallback);
            }
        }
    }
    //! Process sequences from a fasta representation, directly creating seed ins seedMap
    void readSequences(FastaRepresentation const & fasta,
                       std::function<void(TwoBitKmer<TwoBitSeedDataType>, KmerOccurrence, size_t)> const & seedInsertCallback,
                       ParallelProgressBar & pb) {
        std::unique_lock<std::mutex> outputLock(mutexOutput_, std::defer_lock); // dummy for extractSeeds()
        auto stepSize = pb.unprotectedProgressBar().step();
        size_t processedSequences = 0;
        for (auto&& elem : fasta.headerToSequence()) {
            auto& fastaSequence = elem.second;
            auto& sequence = fastaSequence.sequence();
            auto& sequenceName = fastaSequence.sequenceName();
            auto& genomeName = fastaSequence.genomeName();
            auto sequenceID = idMap_->querySequenceIDConst(sequenceName, genomeName);
            auto genomeID = idMap_->queryGenomeIDConst(genomeName);
            extractSeeds(sequence,
                         sequenceID,
                         genomeID,
                         sequenceName,
                         seedInsertCallback, // directly insert
                         outputLock);
            ++processedSequences;
            // update progress bar occasionally
            if (processedSequences >= stepSize) {
                //pb += processedSequences;
                pb.increase(processedSequences);
                processedSequences = 0;
            }
        }
        // final progress bar update if neccessary
        if (processedSequences > 0) {
            pb.increase(processedSequences);
        }
    }
    //! Process sequences from a fasta representation in parallel
    void readSequenceParallel(std::function<void(TwoBitKmer<TwoBitSeedDataType>, KmerOccurrence, size_t)> const & seedInsertCallback,
                              ParallelProgressBar & pb,
                              typename FastaRepresentation::FastaMapType::const_iterator sequencesIt,
                              typename FastaRepresentation::FastaMapType::const_iterator sequencesEnd) {
        // create locks but don't lock yet
        std::unique_lock<std::mutex> memberAccessLock(mutexMemberAccess_, std::defer_lock);
        std::unique_lock<std::mutex> outputLock(mutexOutput_, std::defer_lock);
        std::vector< std::tuple<TwoBitKmer<TwoBitSeedDataType>, KmerOccurrence, size_t> > occurrences;
        auto bufferOccurrence = [&occurrences](TwoBitKmer<TwoBitSeedDataType> seed, KmerOccurrence occ, size_t i) {
            occurrences.emplace_back(seed, occ, i);
        };

        for (; sequencesIt != sequencesEnd; ++sequencesIt) {
            auto& fastaSequence = (*sequencesIt).second;
            auto& sequence = fastaSequence.sequence();
            auto& sequenceName = fastaSequence.sequenceName();
            auto& genomeName = fastaSequence.genomeName();
            auto sequenceID = idMap_->querySequenceIDConst(sequenceName, genomeName);
            auto genomeID = idMap_->queryGenomeIDConst(genomeName);
            extractSeeds(sequence,
                         sequenceID,
                         genomeID,
                         sequenceName,
                         bufferOccurrence,
                         outputLock);
            pb.increase();
        }

        // create seeds
        memberAccessLock.lock();
        for (auto&& tuple : occurrences) {
            seedInsertCallback(std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple));
        }
        memberAccessLock.unlock();

        return;
    }
    //! Helper factory function to check if input data is valid
    void validateInputFiles(FastaCollection const & fastaCollection) const {
        auto genome1Exists = false;
        auto genome2Exists = false;
        std::vector<std::string> inputGenomes;
        for (auto&& elem : fastaCollection.collection()) {
            auto& file = elem.first;
            std::filesystem::path path(file);
            path.replace_extension("");
            inputGenomes.emplace_back(path.filename());
            if (path.filename() == config_->genome1()) { genome1Exists = true; }
            if (path.filename() == config_->genome2()) { genome2Exists = true; }
        }

        if (!(genome1Exists && genome2Exists)) {
            std::cerr << std::endl << "Input genomes: " << inputGenomes << std::endl;
            throw std::runtime_error("[ERROR] -- ExactSeeds -- Genomes not present in input files");
        }
    }

    //! Global configuration
    std::shared_ptr<Configuration const> config_;
    //! Identifier Mapping
    std::shared_ptr<IdentifierMapping const> idMap_;
    //! Masks used for seed creation
    std::shared_ptr<SpacedSeedMaskCollection const> maskCollection_;
    //! Lock for member access of this class
    std::mutex mutexMemberAccess_;
    //! Lock for stdout
    std::mutex mutexOutput_;
    //! For thinning
    std::hash<std::string> stringhasher_;
};

#endif // EXTRACTSEEDS_H
