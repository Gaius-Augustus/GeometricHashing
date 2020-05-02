#include "ExtractSeeds.h"

namespace fs = std::experimental::filesystem;

template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
ExtractSeeds<TwoBitKmerDataType,
             TwoBitSeedDataType>::ExtractSeeds(std::shared_ptr<FastaCollection const> fastaCollection,
                                               std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                                              TwoBitSeedDataType>> seedMap,
                                               int nThreads, bool noProgressbar, bool skipMatchCreation)
    : fastaCollection_{fastaCollection}, mutexMemberAccess_{}, mutexOutput_{},
      noProgressbar_{noProgressbar}, nThreads_{nThreads}, seedMap_{seedMap},
      skipMatchCreation_{skipMatchCreation} {
    // search for genome0 and 1 in inputFiles and abort if not found
    validateInputFiles();
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void ExtractSeeds<TwoBitKmerDataType,
                  TwoBitSeedDataType>::extractAllVsAll() {
    // read files one by one and fill seedMap_
    for (auto&& elem : fastaCollection_->collection()) {
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

    if (skipMatchCreation_) {
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
    std::vector<std::pair<TwoBitKmer<TwoBitKmerDataType>, KmerOccurrence>> occurrences;
    auto span = seedMap_->span();
    auto const idMap = seedMap_->idMapCopy();
    size_t processedSequences = 0;
    // vvv -- attention: uninitialized
    uint16_t genomeID;
    std::string sequence;
    uint32_t sequenceID;
    std::string sequenceName;

    for (; sequencesIt != sequencesEnd; ++sequencesIt) {
        auto& fastaSequence = (*sequencesIt).second;
        genomeID = idMap.queryGenomeIDConst(genomeName);
        sequence = fastaSequence.sequence();
        sequenceName = fastaSequence.sequenceName();
        sequenceID = idMap.querySequenceID(sequenceName);

        extractSeeds(sequence,
                     sequenceID,
                     sequenceName,
                     genomeID,
                     occurrences,
                     span,
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

    // create matches
    memberAccessLock.lock();
    for (auto&& pair : occurrences) {
        seedMapGlobal->createSeed(pair.first, pair.second);
    }
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
                  TwoBitSeedDataType>::extractSeeds(std::string const & sequence,
                                                    uint32_t sequenceID,
                                                    std::string const & sequenceName,
                                                    uint16_t genomeID,
                                                    std::vector<std::pair<TwoBitKmer<TwoBitKmerDataType>, KmerOccurrence>> & occurrences,
                                                    size_t span,
                                                    std::unique_lock<std::mutex> & outputLock) {
    if (sequence.size() < span) {
        outputLock.lock();
        std::cout << "[WARNING] -- ExactKmerMap -- Sequence " << sequenceName
                  << " is shorter than l (" << span << ")" << std::endl;
        outputLock.unlock();
    } else {
        auto validFlag = std::make_shared<bool>();
        for (size_t p = 0; p < (sequence.size() - span + 1); ++p) {
            auto kmer = sequence.substr(p, span);
            TwoBitKmer<TwoBitKmerDataType> twoBitKmer(kmer, validFlag);
            if (*validFlag) {
                auto occurrence = KmerOccurrence(genomeID, sequenceID, p, false, kmer);
                occurrences.emplace_back(twoBitKmer, occurrence);
            }
        }
    }
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void ExtractSeeds<TwoBitKmerDataType,
                  TwoBitSeedDataType>::validateInputFiles() {
    auto genome0Exists = false;
    auto genome1Exists = false;
    std::vector<std::string> inputGenomes;
    for (auto&& elem : fastaCollection_->collection()) {
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



template class ExtractSeeds<TwoBitKmerDataLong, TwoBitKmerDataLong>;
template class ExtractSeeds<TwoBitKmerDataLong, TwoBitKmerDataMedium>;
template class ExtractSeeds<TwoBitKmerDataLong, TwoBitKmerDataShort>;
template class ExtractSeeds<TwoBitKmerDataMedium, TwoBitKmerDataMedium>;
template class ExtractSeeds<TwoBitKmerDataMedium, TwoBitKmerDataShort>;
template class ExtractSeeds<TwoBitKmerDataShort, TwoBitKmerDataShort>;



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void ExtractSeedsFast<TwoBitKmerDataType,
                      TwoBitSeedDataType>::extractPairwise() {
    // all genome pairings and inside each pairing, all sequence pairs
    SequencePairVector sequencePairs;
    auto it1 = this->fastaCollection_->collection().begin();
    auto end = this->fastaCollection_->collection().end();
    while (it1 != end) {
        auto it2 = std::next(it1, 1);
        while (it2 != end) {
            for (auto&& faSeq1 : it1->second.headerToSequence()) {
                for (auto&& faSeq2 : it2->second.headerToSequence()) {
                    sequencePairs.emplace_back(SequencePair{faSeq1.second, faSeq2.second});
                }
            }
            ++it2;
        }
        ++it1;
    }
    std::cout << "Created " << sequencePairs.size() << " sequence pairs, processing..." << std::endl;

    // parallel on chunks of pairs, execute pipeline
    Timestep ts("Extracting seeds, creating matches and running filter on "
                + std::to_string(this->nThreads_) + " threads");
    ProgressBar pb(sequencePairs.size(), this->noProgressbar_);

    // spawn threads
    std::vector<std::thread> threads;
    for (auto i = 0; i < this->nThreads_; ++i) {  // using FastaMapType = std::unordered_map<std::string, FastaSequence>;
        auto chunkIt = getContainerChunkBegin<SequencePairVector>(sequencePairs, i, this->nThreads_);
        auto chunkEnd = getContainerChunkEnd<SequencePairVector>(sequencePairs, i, this->nThreads_);
        if (chunkIt != chunkEnd) {
            threads.push_back(std::thread(&ExtractSeedsFast<TwoBitKmerDataType,
                                                            TwoBitSeedDataType>::runPipeline, this,
                                          chunkIt, chunkEnd, std::ref(pb)));
        }
    }

    // wait until all threads are finished
    std::for_each(threads.begin(), threads.end(), [](std::thread & t) { t.join(); });
    std::cout << std::endl;

    ts.endAndPrint(Timestep::minutes);
    std::cout << std::endl;
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
void ExtractSeedsFast<TwoBitKmerDataType,
                      TwoBitSeedDataType>::runPipeline(typename SequencePairVector::const_iterator pairsIt,
                                                       typename SequencePairVector::const_iterator pairsEnd,
                                                       ProgressBar & pb) {
    // create locks but don't lock yet
    std::unique_lock<std::mutex> outputLock(mutexOutfile_, std::defer_lock);
    std::unique_lock<std::mutex> stdoutLock(mutexStdout_, std::defer_lock);
    auto stepSize = pb.step();
    auto seedMap = this->seedMap_->localCopy();
    auto osbuf = std::stringstream();
    std::vector<std::pair<TwoBitKmer<TwoBitKmerDataType>, KmerOccurrence>> occurrences;
    auto span = this->seedMap_->span();

    // perform entire pipeline for each pair
    size_t pairsProcessed = 0;
    size_t batchProgress = 0;
    size_t totalProcessed = 0;
    for (; pairsIt != pairsEnd; ++pairsIt) {
        // extract seeds
        for (auto&& fastaSequenceRef : *pairsIt) {
            auto& fastaSequence = fastaSequenceRef.get();
            auto const & genomeName = fastaSequence.genomeName();
            auto genomeID = seedMap->idMap()->queryGenomeIDConst(genomeName);
            auto const & sequence = fastaSequence.sequence();
            auto const & sequenceName = fastaSequence.sequenceName();
            auto sequenceID = seedMap->idMap()->querySequenceID(sequenceName);
            this->extractSeeds(sequence,
                               sequenceID,
                               sequenceName,
                               genomeID,
                               occurrences,
                               span,
                               stdoutLock);
        }
        ++batchProgress; ++pairsProcessed; ++totalProcessed;

        // (last) batch full OR (if batches AND batch size reached)
        if (std::next(pairsIt, 1) == pairsEnd
                || (config_->fastBatchsize() > 0
                    && batchProgress >= config_->fastBatchsize())) {
            stdoutLock.lock();
            std::cout << "Thread " << std::this_thread::get_id() << " reached batch size, starting pipeline" << std::endl;
            stdoutLock.unlock();
            // create seedMap_
            for (auto&& pair : occurrences) {
                seedMap->createSeed(pair.first, pair.second);
            }
            occurrences.clear();
            // create matches if appropriate
            if (!config_->performGeometricHashing()) {
                seedMap->createMatches();
                stdoutLock.lock();
                std::cout << "Thread " << std::this_thread::get_id() << " created matches" << std::endl;
                stdoutLock.unlock();
            }

            // perform filtering if appropriate
            if (config_->performGeometricHashing()) {
                // this includes M4 filtering inside Cubes
                geometricHashing<TwoBitKmerDataType, TwoBitSeedDataType>(seedMap->idMap(),
                                                                         seedMap,
                                                                         this->fastaCollection_,
                                                                         config_,
                                                                         false); // no cube output
                stdoutLock.lock();
                std::cout << "Thread " << std::this_thread::get_id() << " finished GH" << std::endl;
                stdoutLock.unlock();
            } else if (config_->performDiagonalFiltering()) {
                seedMap->applyDiagonalMatchesFilter(std::make_shared<std::array<size_t, 4>>());
                stdoutLock.lock();
                std::cout << "Thread " << std::this_thread::get_id() << " finished M4" << std::endl;
                stdoutLock.unlock();
            }

            // writing matches to outfile (or buffer if batchsize 1)
            if (config_->fastBatchsize() == 1) {
                seedMap->output(osbuf);
            } else {
                outputLock.lock();
                seedMap->output(*os_);
                outputLock.unlock();
            }

            // update progress bar occasionally
            if (!this->noProgressbar_ && (pairsProcessed >= stepSize)) {
                stdoutLock.lock();
                pb += pairsProcessed;
                stdoutLock.unlock();
                pairsProcessed = 0;
            }

            // cleanup
            seedMap->clearMatches();
            seedMap->clearSeedMap();

            batchProgress = 0;

            stdoutLock.lock();
            std::cout << "Thread " << std::this_thread::get_id() << " finished pipeline, processed " << totalProcessed << " pairs" << std::endl;
            stdoutLock.unlock();
        }
    }

    // final progress bar update if neccessary
    if (!this->noProgressbar_ && (pairsProcessed > 0)) {
        stdoutLock.lock();
        pb += pairsProcessed;
        stdoutLock.unlock();
    }

    // actual write to file if batchsize 1
    if (config_->fastBatchsize() == 1) {
        outputLock.lock();
        *os_ << osbuf.rdbuf() ;
        outputLock.unlock();
    }
}



template class ExtractSeedsFast<TwoBitKmerDataLong, TwoBitKmerDataLong>;
template class ExtractSeedsFast<TwoBitKmerDataLong, TwoBitKmerDataMedium>;
template class ExtractSeedsFast<TwoBitKmerDataLong, TwoBitKmerDataShort>;
template class ExtractSeedsFast<TwoBitKmerDataMedium, TwoBitKmerDataMedium>;
template class ExtractSeedsFast<TwoBitKmerDataMedium, TwoBitKmerDataShort>;
template class ExtractSeedsFast<TwoBitKmerDataShort, TwoBitKmerDataShort>;
