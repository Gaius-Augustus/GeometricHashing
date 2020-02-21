#include "ExtractSeeds.h"

namespace fs = std::experimental::filesystem;

template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
ExtractSeeds<TwoBitKmerDataType,
             TwoBitSeedDataType>::ExtractSeeds(std::shared_ptr<FastaCollection const> fastaCollection,
                                               std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                                              TwoBitSeedDataType>> seedMap,
                                               int nThreads, bool noProgressbar, bool skipMatchCreation)
    : mutexMemberAccess_{}, mutexOutput_{},
      noProgressbar_{noProgressbar}, nThreads_{nThreads},
      seedMap_{seedMap} {
    // search for genome0 and 1 in inputFiles and abort if not found
    validateInputFiles(fastaCollection);
    // run the appropriate extraction
    extractAllVsAll(fastaCollection, skipMatchCreation);
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void ExtractSeeds<TwoBitKmerDataType,
                  TwoBitSeedDataType>::extractAllVsAll(std::shared_ptr<FastaCollection const> fastaCollection, bool skipMatchCreation) {
    // read files one by one and fill seedMap_
    for (auto&& elem : fastaCollection->collection()) {
        auto& fasta = elem.second;

        Timestep ts("Extracting k-mers from " + fasta.genome() + " on "
                    + std::to_string(nThreads_) + " threads");
        ProgressBar pb(fasta.numSequences(), noProgressbar_);

        // spawn threads
        std::vector<std::thread> threads;

        for (auto i = 0; i < nThreads_; ++i) {  // using FastaMapType = std::unordered_map<std::string, FastaSequence>;
            auto chunkIt = getContainerChunkBegin<FastaRepresentation::FastaMapType>(fasta.headerToSequence(),
                                                                                     i, nThreads_);
            auto chunkEnd = getContainerChunkEnd<FastaRepresentation::FastaMapType>(fasta.headerToSequence(),
                                                                                    i, nThreads_);
            if (chunkIt != chunkEnd) {
                // https://thispointer.com/c11-start-thread-by-member-function-with-arguments/
                threads.push_back(std::thread(&ExtractSeeds<TwoBitKmerDataType,
                                                            TwoBitSeedDataType>::readSequence, this,
                                              seedMap_,
                                              chunkIt, chunkEnd, fasta.genome(), std::ref(pb)));
            }
        }

        // wait until all threads are finished
        std::for_each(threads.begin(), threads.end(), [](std::thread & t) { t.join(); });
        std::cout << std::endl;

        ts.endAndPrint(Timestep::seconds);
        std::cout << std::endl;
    }

    seedMap_->printStatistics();

    if (skipMatchCreation) {
        std::cout << "[INFO] -- ExtractSeeds::extractAllVsAll -- Skipping match creation" << std::endl;
    } else {
        Timestep ts("Creating matches");
        seedMap_->createMatchesGeneral();
        ts.endAndPrint(Timestep::minutes);

        std::cout << "Created " << seedMap_->numMatches() << " matches" << std::endl;
    }
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void ExtractSeeds<TwoBitKmerDataType,
                  TwoBitSeedDataType>::extractSeeds(std::string const & sequence,
                                                    uint32_t sequenceID,
                                                    std::string const & sequenceName,
                                                    uint16_t genomeID,
                                                    std::string const & genomeName,
                                                    std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> & seedMapLocal,
                                                    std::unique_lock<std::mutex> & outputLock) {
    if (sequence.size() < seedMapLocal->span()) {
        outputLock.lock();
        std::cout << "[WARNING] -- ExactKmerMap -- Sequence " << sequenceName
                  << " is shorter than k (" << seedMapLocal->span() << ")" << std::endl;
        outputLock.unlock();
    } else {
        seedMapLocal->updateSequenceLengthStatistics(genomeName, sequence.size());
        auto validFlag = std::make_shared<bool>();

        for (size_t p = 0; p < (sequence.size() - seedMapLocal->span() + 1); ++p) {
            auto kmer = sequence.substr(p, seedMapLocal->span());

            TwoBitKmer<TwoBitKmerDataType> twoBitKmer(kmer, validFlag);
            if (*validFlag) {
                auto occurrence = KmerOccurrence(genomeID, sequenceID, p, false, kmer);
                seedMapLocal->createSeed(twoBitKmer, occurrence); // handle filtering there
            } else {
                seedMapLocal->updateInvalidSeedStatistics(genomeName);
            }
        }
    }
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void ExtractSeeds<TwoBitKmerDataType,
                  TwoBitSeedDataType>::readSequence(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> seedMapGlobal,
                                                    typename FastaRepresentation::FastaMapType::const_iterator sequencesIt,
                                                    typename FastaRepresentation::FastaMapType::const_iterator sequencesEnd,
                                                    std::string genomeName,
                                                    ProgressBar & pb) {
    // create locks but don't lock yet
    std::unique_lock<std::mutex> memberAccessLock(mutexMemberAccess_, std::defer_lock);
    std::unique_lock<std::mutex> outputLock(mutexOutput_, std::defer_lock);
    auto stepSize = pb.step();
    // store data local first (avoid excessive blocking among threads)
    memberAccessLock.lock();
    auto seedMapLocal = seedMapGlobal->localCopy();
    memberAccessLock.unlock();

    size_t processedSequences = 0;

    // vvv -- attention: uninitialized
    uint16_t genomeID;
    std::string sequence;
    uint32_t sequenceID;
    std::string sequenceName;

    for (; sequencesIt != sequencesEnd; ++sequencesIt) {
        auto& fastaSequence = (*sequencesIt).second;
        genomeID = seedMap_->idMap()->queryGenomeIDConst(genomeName);
        sequence = fastaSequence.sequence();
        sequenceName = fastaSequence.sequenceName();
        sequenceID = seedMap_->idMap()->querySequenceID(sequenceName);

        extractSeeds(sequence,
                     sequenceID,
                     sequenceName,
                     genomeID,
                     genomeName,
                     seedMapLocal,
                     outputLock);
        ++processedSequences;

        // update progress bar occasionally
        if (!noProgressbar_ && (processedSequences >= stepSize)) {
            outputLock.lock();
            pb += processedSequences;
            outputLock.unlock();
            processedSequences = 0;
        }
    }

    // merge local data back
    memberAccessLock.lock();
    seedMapGlobal->merge(seedMapLocal);
    memberAccessLock.unlock();

    // final progress bar update if neccessary
    if (!noProgressbar_ && processedSequences > 0) {
        outputLock.lock();
        pb += processedSequences;
        outputLock.unlock();
    }

    return;
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void ExtractSeeds<TwoBitKmerDataType,
                  TwoBitSeedDataType>::validateInputFiles(std::shared_ptr<FastaCollection const> fastaCollection) {
    auto genome0Exists = false;
    auto genome1Exists = false;
    std::vector<std::string> inputGenomes;
    for (auto&& elem : fastaCollection->collection()) {
        auto& file = elem.first;
        fs::path path(file);
        path.replace_extension("");
        inputGenomes.emplace_back(path.filename());
        if (path.filename() == genome0()) { genome0Exists = true; }
        if (path.filename() == genome1()) { genome1Exists = true; }
    }

    if (!(genome0Exists && genome1Exists)) {
        std::cout << std::endl << "Input genomes: " << inputGenomes << std::endl;
        throw std::runtime_error("[ERROR] -- ExactSeeds -- Genomes not present in input files");
    }
}



// Normally not possible to separate definitions of a template class from declarations
// However, workaround: have compiler create instantiations of both possible TwoBitKmerTypes
// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
template class ExtractSeeds<TwoBitKmerDataLong, TwoBitKmerDataLong>;
template class ExtractSeeds<TwoBitKmerDataLong, TwoBitKmerDataMedium>;
template class ExtractSeeds<TwoBitKmerDataLong, TwoBitKmerDataShort>;
template class ExtractSeeds<TwoBitKmerDataMedium, TwoBitKmerDataMedium>;
template class ExtractSeeds<TwoBitKmerDataMedium, TwoBitKmerDataShort>;
template class ExtractSeeds<TwoBitKmerDataShort, TwoBitKmerDataShort>;
