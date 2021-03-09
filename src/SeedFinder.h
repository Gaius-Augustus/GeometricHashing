#ifndef SEEDFINDER_H
#define SEEDFINDER_H

#include <climits>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>

#include "Configuration.h"
#include "Cubeset.h"
#include "DiagonalMatchesFilter.h"
#include "ExtractSeeds.h"
#include "FastaCollection.h"
#include "IdentifierMapping.h"
#include "Linkset.h"
#include "Output.h"
#include "ParallelizationUtils.h"
#include "ParallelProgressBarHandler.h"
#include "SeedMap.h"
#include "SpacedSeedMaskCollection.h"
#include "TwoBitKmer.h"

// Create SeedMap or SeedKmerMap based on the template parameters

//! Pass information of verbosity level and whether it is safe/desired to run something in parallel
struct ParallelVerboseInfo {
    ParallelVerboseInfo(bool allowParallelExecution, bool zeroOutput)
        : allowParallelExecution{allowParallelExecution}, zeroOutput{zeroOutput} {}
    bool allowParallelExecution;
    bool zeroOutput;
};




// Some code to tell SeedMap and SeedKmerMap apart
// default: whatever you pass is not a SeedKmerMap
inline constexpr auto isSeedKmerMapImpl(...) { return std::false_type{}; }

// if M has IsSeedKmerMap type, say it is a SeedKmerMap
template<typename M, typename = typename M::IsSeedKmerMap>
constexpr auto isSeedKmerMapImpl(M const &) { return std::true_type{}; }

// Call this to tell SeedMap and SeedKmerMap apart (returns true_type or false_type)
template<typename SeedMapType>
constexpr auto isSeedKmerMap(SeedMapType const & m) { return isSeedKmerMapImpl(m); }

// Now, factory function templates depending on above checks
template<typename SeedMapType>
auto createSeedMapImpl(std::shared_ptr<SeedMapType> seedMap,
                       std::shared_ptr<FastaCollection const> fastaCollection,
                       ParallelVerboseInfo const & pinf,
                       std::false_type) { // for SeedMap
    Timestep tsDirectExtract("Extracting seeds from input files", pinf.zeroOutput);
    seedMap->extractSeeds(fastaCollection, pinf.allowParallelExecution);
    seedMap->cleanupReferenceSeeds(pinf.allowParallelExecution);
    tsDirectExtract.endAndPrint();
    if (!pinf.zeroOutput) { seedMap->printStatistics(); }
    return seedMap;
}

template<typename SeedMapType>
auto createSeedMapImpl(std::shared_ptr<SeedMapType> seedKmerMap,
                       std::shared_ptr<FastaCollection const> fastaCollection,
                       ParallelVerboseInfo const & pinf,
                       std::true_type) { // for SeedKmerMap
    (void)seedKmerMap; (void)fastaCollection; (void)pinf;
    throw std::runtime_error("[ERROR] -- createSeedMapImpl -- SeedKmerMap not available in this version");
}



//! Create partial FastaCollection from a complete FastaCollection and a vector of sequence IDs
/*! The sequence IDs must denote a subset of sequences present in \c completeFastaCollection */
inline std::shared_ptr<FastaCollection const> getPartialFastaCollection(std::vector<size_t> const & sequenceIDs,
                                                                        FastaCollection const & completeFastaCollection,
                                                                        IdentifierMapping const & idMap) {
    auto fastaCollection = std::make_shared<FastaCollection>();
    tsl::hopscotch_map<std::string, FastaRepresentation> genomeToFasta;
    for (auto sid : sequenceIDs) {
        auto& seqTuple = idMap.querySequenceTuple(sid);
        auto& genome = idMap.queryGenomeName(seqTuple.gid);
        if (genomeToFasta.find(genome) == genomeToFasta.end()) {
            genomeToFasta.insert({genome, FastaRepresentation{FastaGenomeName{genome}}});
        }
        auto& faSeq = completeFastaCollection.fastaSequence(sid, idMap);
        genomeToFasta.at(genome).addSequence(faSeq.sequenceName(), faSeq.sequence());
    }
    for (auto&& elem : genomeToFasta) {
        fastaCollection->emplace(elem.first, elem.second);
    }
    return fastaCollection;
}



//! Implementing different calling scenarios of the basic pipeline
/*! The basic pipeline consists of creeating links from a SeedMap,
 * possibly filtering those links via M4 or GH, grouping that links
 * and writing them to the output file.
 *
 * Calling Sceanrios are all-vs-all sequences at once, either as
 *   single run or multiple runs on different SeedMaps in parallel
 *   and 1-vs-all (one reference genome sequence at a time against
 *   all non-ref-sequences) either contiguously or in parallel,
 *   each either with or without M4/GH filtering.
 *   A special scenario is the two-step GH */
template<typename SeedMapType, typename LinksetType>
class BasicPipeline {
public:
    static struct AllVsAll{} allVsAll;
    static struct OneVsAll{} oneVsAll;
    static struct PreFilter{} preFilter;

    BasicPipeline(std::shared_ptr<FastaCollection const> fastaCollection,
                  std::shared_ptr<Output> output,
                  std::shared_ptr<IdentifierMapping const> idMap,
                  std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> seqLens,
                  std::shared_ptr<Configuration const> config)
        : config_{config}, fastaCollection_{fastaCollection},
          idMap_{idMap}, mutexOutput_{}, output_{output},
          seqLens_{seqLens} {}

    void run(AllVsAll, ParallelVerboseInfo const & pinf) {
        Timestep tsSeedMap("~~~ Create Seed Map (all-vs-all) ~~~", pinf.zeroOutput);
        auto seedMap = createSeedMap(fastaCollection_, pinf); // if run() allowed to start parallel execution, run seed extraction in parallel
        tsSeedMap.endAndPrint();
        Timestep ts("~~~ Create Linkset and Output Matches (all-vs-all) ~~~", pinf.zeroOutput); // if silent, timestep in silent mode
        MM::MemoryMonitor mm;
        if (!pinf.zeroOutput) { std::cout << "Memory usage before link creation" << std::endl << mm << std::endl; }
        Timestep tsLinkset{"Creating Links", pinf.zeroOutput};

        auto linkset = std::make_shared<LinksetType>(config_, seedMap->idMap(), pinf.allowParallelExecution);
        linkset->createLinks(*seedMap, pinf.zeroOutput);

        if (!pinf.zeroOutput) { std::cout << "Memory usage after link creation" << std::endl << mm << std::endl; }

        seedMap->clear();    // save memory

        if (!pinf.zeroOutput) {
            std::cout << "[INFO] -- Created " << linkset->size() << " links" << std::endl;
            std::cout << "Memory usage after clearing seedMap" << std::endl << mm << std::endl;
        }
        tsLinkset.endAndPrint();

        if (config_->performDiagonalFiltering()) {
            Timestep tsDiag("Apply Diagonal Filter", pinf.zeroOutput);
            linkset->applyDiagonalMatchesFilter();
            tsDiag.endAndPrint();
        } else if (config_->performGeometricHashing()) {
            Timestep tsGH("Run Geometric Hashing", pinf.zeroOutput);
            geometricHashing(linkset, seqLens_, pinf);
            tsGH.endAndPrint();
        } else {
            Timestep tsGroup("Grouping overlapping seeds", pinf.zeroOutput);
            linkset->groupOverlappingLinks();
            tsGroup.endAndPrint();
        }

        if (!pinf.zeroOutput) { linkset->printStatistics(std::cout); }
        Timestep tsOutput("Writing matches to output", pinf.zeroOutput);

        output_->outputLinkset(*linkset, pinf.zeroOutput);

        tsOutput.endAndPrint();
        ts.endAndPrint();
    }
    void run(OneVsAll, ParallelVerboseInfo const & pinf) {
        Timestep tsSeedMap("~~~ Create Seed Map (one-vs-all) ~~~", pinf.zeroOutput);
        auto seedMap = createSeedMap(fastaCollection_, pinf);
        tsSeedMap.endAndPrint();
        auto fun = [this,
                    &seedMap](ParallelVerboseInfo lambdaPinf,
                              typename std::vector<size_t>::const_iterator it,
                              typename std::vector<size_t>::const_iterator end) {
            auto linkset = std::make_shared<LinksetType>(config_, seedMap->idMap(), lambdaPinf.allowParallelExecution);
            for (; it != end; ++it) {
                linkset->createLinks(*seedMap, *it);
                if (config_->performDiagonalFiltering()) {
                    linkset->applyDiagonalMatchesFilter();
                } else if (config_->performGeometricHashing()) {
                    geometricHashing(linkset, seqLens_, lambdaPinf);
                } else {
                    linkset->groupOverlappingLinks();
                }
                if (output_->tryOutputLinkset(*linkset, lambdaPinf.zeroOutput)) { linkset->clear(); }
            }
            if (linkset->size()) {
              output_->outputLinkset(*linkset, lambdaPinf.zeroOutput);
              linkset->clear();
            }
        };
        std::vector<size_t> sids;
        for (auto&& elem : seedMap->referenceSeedMap().referenceSeedMap()) { sids.emplace_back(elem.first); }
        if (pinf.allowParallelExecution && config_->nThreads() > 1) {
            executeParallel(sids, config_->nThreads(), fun, ParallelVerboseInfo{false, false}); // no extra parallel execution and no output
        } else {
            executeParallel(sids, 1, fun, pinf); // no extra thread spawned, forward pinf
        }
    }
    void run(AllVsAll, PreFilter, ParallelVerboseInfo const & pinf) {
        Timestep tsSeedMap("~~~ Create Pre-Seed Map (all-vs-all) ~~~", pinf.zeroOutput);
        auto seedMap = createSeedMap(fastaCollection_, pinf, true); // get pre-GH seedMap
        tsSeedMap.endAndPrint();
        Timestep tsPre("~~~ Find Relevant Sequence Tuples/Cubes (all-vs-all) ~~~", pinf.zeroOutput);
        MM::MemoryMonitor mm;
        if (!pinf.zeroOutput) { std::cout << "Memory usage before pre-filter" << std::endl << mm << std::endl; }
        PrefilterCubeset preCubeset{*seedMap, idMap_, seqLens_, config_, pinf.zeroOutput};
        if (!pinf.zeroOutput) { std::cout << "Memory usage after pre-filter" << std::endl << mm << std::endl; }
        seedMap->clear(); // save memory
        tsPre.endAndPrint();
        if (!pinf.zeroOutput) { std::cout << "[INFO] -- " << preCubeset.relevantCubeSet().size()
                                          << " relevant seqtuples/cubes form " << preCubeset.sequenceCluster().size()
                                          << " clusters" << std::endl; }
        if (config_->cubeOutput() > 0) {
            if (!pinf.zeroOutput) { std::cout << "[INFO] -- Writing relevant cubes to disk" << std::endl; }
            output_->outputRelevantCubes(preCubeset); // output locks all its public functions, thread safe
        }

        Timestep tsPostGH("~~~ Run Seed Finding on Sequence Tuples/Relevant Cubes ~~~", pinf.zeroOutput);
        auto sequenceCluster = preCubeset.sequenceCluster();
        ParallelProgressBar pb{sequenceCluster.size(), pinf.zeroOutput || config_->verbose() < 2};
        // for each cluster, run post-GH
        auto runPostGH = [this, &pb](ParallelVerboseInfo lambdaPinf,
                                     std::vector<std::shared_ptr<SequenceCluster>>::const_iterator it,
                                     std::vector<std::shared_ptr<SequenceCluster>>::const_iterator end) {
            Timestep tsBatch("Running GH for batch", lambdaPinf.zeroOutput);
            for (; it != end; ++it) {
                auto& cluster = **it;
                std::vector<size_t> sidv{cluster.sids.begin(), cluster.sids.end()};
                auto partialFastaCollection = getPartialFastaCollection(sidv, *fastaCollection_, *idMap_);
                //  -> idMapping needs to be consistent, only recreate seqLens for Cube scoring
                auto partialSequenceLengths = std::make_shared<tsl::hopscotch_map<size_t, size_t>>();
                partialFastaCollection->fillSequenceLengths(*partialSequenceLengths, *idMap_);
                // get SeedMap (overwrite pre-seedMap)
                Timestep tsPostSeedMap("Extracting Seeds for Second Run", lambdaPinf.zeroOutput);
                auto seedMap = createSeedMap(partialFastaCollection, lambdaPinf, false); // use regular seed masks
                tsPostSeedMap.endAndPrint();
                // create linkset with only relevant links
                Timestep tsLinkset("Create Relevant Links", lambdaPinf.zeroOutput);
                auto linkset = std::make_shared<LinksetType>(config_, idMap_, lambdaPinf.allowParallelExecution);
                //linkset->createLinks(*seedMap, preCubeset.relevantCubeSet());
                linkset->createLinks(*seedMap, cluster.cubes, lambdaPinf.zeroOutput);
                seedMap->clear(); // save memory
                tsLinkset.endAndPrint();
                if (config_->performGeometricHashing()) {
                    Timestep tsGH("Run Geometric Hashing", lambdaPinf.zeroOutput);
                    geometricHashing(linkset, partialSequenceLengths, lambdaPinf);
                    tsGH.endAndPrint();
                } else {
                    Timestep tsGroup("Grouping overlapping seeds", lambdaPinf.zeroOutput);
                    linkset->groupOverlappingLinks();
                    tsGroup.endAndPrint();
                }
                if (!lambdaPinf.zeroOutput) { linkset->printStatistics(std::cout); }
                Timestep tsOutput("Writing matches to output", lambdaPinf.zeroOutput);
                output_->outputLinkset(*linkset, lambdaPinf.zeroOutput);
                tsOutput.endAndPrint();
                pb.increase();
            }
            tsBatch.endAndPrint();
        };
        if ((!config_->postSequential()) && pinf.allowParallelExecution && config_->nThreads() > 1) {
            executeParallel(sequenceCluster, config_->nThreads(), runPostGH, ParallelVerboseInfo{false, true}); // no extra parallel execution and no output
        } else {
            executeParallel(sequenceCluster, 1, runPostGH, pinf); // no extra thread spawned, forward pinf
        }
        tsPostGH.endAndPrint();
    }

private:
    std::shared_ptr<SeedMapType> createSeedMap(std::shared_ptr<FastaCollection const> fastaCollection,
                                               ParallelVerboseInfo const & pinf,
                                               bool preFilter = false) const {
        auto seedMap = (preFilter)
                ? std::make_shared<SeedMapType>(config_, idMap_, PrefilterSeedMapTag{})
                : std::make_shared<SeedMapType>(config_, idMap_);
        return createSeedMapImpl(seedMap,
                                 fastaCollection,
                                 pinf,
                                 isSeedKmerMap(*seedMap));
    }
    //! Normal GH pipeline
    void geometricHashing(std::shared_ptr<LinksetType> const & linkset,
                          std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> const & seqLens,
                          ParallelVerboseInfo const & pinf) {
        std::unique_lock<std::mutex> outputLock(mutexOutput_, std::defer_lock);
        Timestep tsCreate("Creating Cubeset", pinf.zeroOutput);
        Cubeset cubeset(linkset, seqLens, pinf.allowParallelExecution);
        if (!pinf.zeroOutput) { std::cout << "Created " << cubeset.cubeMap().size() << " Cubes" << std::endl; }
        tsCreate.endAndPrint();
        linkset->clear();   // save memory
        Timestep tsExtract("Extracting Matches from Cubes", pinf.zeroOutput);
        if (config_->cubeOutput() > 0) {
            output_->outputCubes(cubeset); // output locks all its public functions, thread safe
        }
        for (auto&& elem : cubeset.scoreToCube()) {
            auto score = elem.first;
            auto& cubes = elem.second;
            if (score >= static_cast<double>(config_->cubeScoreThreshold())) {
                for (auto&& cubeptr : cubes) {
                    for (auto&& link : cubeset.cubeMap().at(cubeptr)) {
                        linkset->createLinks(link.occurrence(), link.span());
                    }
                    if (config_->hasse()) {
                        if (cubeset.subcubeMap().find(cubeptr) != cubeset.subcubeMap().end()) {
                            for (auto subcube : cubeset.subcubeMap().at(cubeptr)) {
                                for (auto&& link : cubeset.cubeMap().at(subcube)) {
                                    linkset->createLinks(link.occurrence(), link.span());
                                }
                            }
                        }
                    }
                }
            }
        }
        linkset->groupOverlappingLinks();
        if (!pinf.zeroOutput) { std::cout << "Created " << linkset->size() << " matches" << std::endl; }
        tsExtract.endAndPrint();
    }

    std::shared_ptr<Configuration const> const config_;
    std::shared_ptr<FastaCollection const> const fastaCollection_;
    std::shared_ptr<IdentifierMapping const> const idMap_;
    std::mutex mutexOutput_;
    std::shared_ptr<Output> const output_;
    std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> const seqLens_;

};



//! Main class of the program
/*! Runs the seed finding pipeline with additional filtering steps, according to the configuration */
template <typename TwoBitSeedDataType>
class SeedFinder {
public:
    //! c'tor, initializes all neccessary data structures to run the program
    /*! \param config Configuration object
     * \details Remember to fill the idMap with all required genome and sequence IDs
     *  before running seed extraction. */
    SeedFinder(std::shared_ptr<Configuration> config)
        : allvsall_{config->allvsall()},
          config_{config}, mutexOutput_{}, output_{std::make_shared<Output>(config_)},
          parallelSeedMap_{true},
          pipelineA_{!(config_->performDiagonalFiltering() || config_->performGeometricHashing())} {
        std::cout << "[INFO] -- Masks used: " << *(config_->maskCollection()) << std::endl << std::endl;
    }

    //! Run the program. Only call this once.
    void run() {
        Timestep tsRun{"Run"};
        auto batchvsbatch = (config_->batchsize() > 1) ? true : false;
        auto fastaInput = config_->inputFiles().size() ? true : false;
        if (fastaInput) {
            std::cout << "[INFO] -- Running with Fasta files as input" << std::endl;
        } else {
            std::cout << "[INFO] -- Running with Metagraph as input" << std::endl;
        }
        if (allvsall_) {
            std::cout << "[INFO] -- Running all-vs-all mode" << std::endl;
        } else {
            std::cout << "[INFO] -- Running 1-vs-all mode" << std::endl;
        }
        if (config_->performGeometricHashing()) {
            // this includes diagonal filtering inside Cubes if set
            std::cout << "[INFO] -- Running with GeometricHashing Filter" << std::endl;
            if (config_->preMaskCollection()) {
                std::cout << "[INFO] -- Running with Pre-GeometricHashing Step" << std::endl;
            }
        }
        if (config_->performDiagonalFiltering()) {
            std::cout << "[INFO] -- Running with Diagonal Matches Filter" << std::endl;
        }

        // global sequence information
        auto completeIDMap = blankIDMap();
        auto completeSequenceLengths = std::make_shared<tsl::hopscotch_map<size_t, size_t>>();
        auto completeFastaCollection = (fastaInput) ? loadFastas() : std::make_shared<FastaCollection>();
        if (fastaInput) {
            fillIDMapSequenceLengths(*completeFastaCollection, *completeIDMap, *completeSequenceLengths);
        } else {
            throw std::runtime_error("[ERROR] -- SeedFinder::run() -- Metagraph not available in this version");
        }

        // RUN ALL-VS-ALL OR 1-VS-ALL IN BATCH MODE
        if (batchvsbatch) {
            std::cout << "[INFO] -- Running batch mode" << std::endl;
            if (!fastaInput) { throw std::runtime_error("[ERROR] -- SeedFinder::run -- cannot run batchmode with metagraph"); }
            // Run Pipeline for each batch
            auto batches = getFastaBatches(*completeIDMap);
            // run pipeline sequentially
            Timestep tsBatch("Running batches sequentially");
            size_t run = 1; size_t total = batches.size();
            for (auto&& sids : batches) {
                Timestep tsBatchN("Running batch ("+std::to_string(run)+"/"+std::to_string(total)+")");
                auto fastaCollection = getPartialFastaCollection(sids, *completeFastaCollection, *completeIDMap);
                auto idMap = blankIDMap();
                auto sequenceLengths = std::make_shared<tsl::hopscotch_map<size_t, size_t>>();
                fillIDMapSequenceLengths(*fastaCollection, *idMap, *sequenceLengths); // only make sequences and IDs of this batch visible (important for GH scoring)
                callBasicPipeline<SeedMap<TwoBitSeedDataType>>(fastaCollection, idMap, sequenceLengths);
                tsBatchN.endAndPrint();
                ++run;
                tsBatch.endAndPrint();
            }
        }
        // RUN ALL-VS-ALL OR 1-VS-ALL ON COMPLETE DATA
        else {
            if (fastaInput) {
                callBasicPipeline<SeedMap<TwoBitSeedDataType>>(completeFastaCollection, completeIDMap, completeSequenceLengths);
            } else {
                throw std::runtime_error("[ERROR] -- SeedFinder::run() -- SeedKmerMap not available in this version");
            }
        }
        tsRun.endAndPrint();
        output_->addRunInfo("runtime", tsRun.elapsed(Timestep::minutes));
    }
    // Getter functions
    auto config() const { return config_; }
    auto fastaBatches(IdentifierMapping const & idMap) const { return getFastaBatches(idMap); } // for testing

private:
    //! Empty idMap, only genome IDs 0 and 1 are set
    std::shared_ptr<IdentifierMapping> blankIDMap() const {
        auto idMap = std::make_shared<IdentifierMapping>(config_->genome1());
        idMap->queryGenomeID(config_->genome2());
        return idMap;
    }
    template<typename SeedMapType>
    void callBasicPipeline(std::shared_ptr<FastaCollection const> fastaCollection,
                           std::shared_ptr<IdentifierMapping const> idMap,
                           std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> sequenceLengths) {
        if (config_->performGeometricHashing()) {
            auto pipeline = BasicPipeline<SeedMapType,
                                          Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>{fastaCollection, output_, idMap, sequenceLengths, config_};
            if (config_->preMaskCollection()) {
                if (allvsall_) {
                    pipeline.run(pipeline.allVsAll, pipeline.preFilter, ParallelVerboseInfo{true, (config_->verbose() == 0)}); // run all-vs-all, never called in parallel thus allow parallel and output if verbose >= 1
                } else {
                    throw std::runtime_error("[ERROR] -- SeedFinder::callBasicPipeline -- Pre-GH in 1-vs-all mode currently not supported");
                }
            } else {
                if (allvsall_) {
                    pipeline.run(pipeline.allVsAll, ParallelVerboseInfo{true, (config_->verbose() == 0)}); // run all-vs-all, never called in parallel thus allow parallel and output if verbose >= 1
                } else {
                    pipeline.run(pipeline.oneVsAll, ParallelVerboseInfo{true, (config_->verbose() == 0)}); // run 1-vs-all, not called in parallel thus allow parallel and output if verbose >= 1
                }
            }
        } else { // no need for link ptrs, save some memory
            BasicPipeline<SeedMapType,
                          Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>> pipeline{fastaCollection, output_, idMap, sequenceLengths, config_};
            if (config_->preMaskCollection()) {
                if (allvsall_) {
                    pipeline.run(pipeline.allVsAll, pipeline.preFilter, ParallelVerboseInfo{true, (config_->verbose() == 0)}); // run all-vs-all, never called in parallel thus allow parallel and output if verbose >= 1
                } else {
                    throw std::runtime_error("[ERROR] -- SeedFinder::callBasicPipeline -- Pre-GH in 1-vs-all mode currently not supported");
                }
            } else {
                if (allvsall_) {
                    pipeline.run(pipeline.allVsAll, ParallelVerboseInfo{true, (config_->verbose() == 0)}); // run all-vs-all, never called in parallel thus allow parallel and output if verbose >= 1
                } else {
                    pipeline.run(pipeline.oneVsAll, ParallelVerboseInfo{true, (config_->verbose() == 0)}); // run 1-vs-all, not called in parallel thus allow parallel and output if verbose >= 1
                }
            }
        }
    }
    //! Fill idMap and sequenceLengths map from a FastaCollection
    void fillIDMapSequenceLengths(FastaCollection const & fastaCollection,
                                  IdentifierMapping & idMap,
                                  tsl::hopscotch_map<size_t, size_t> & sequenceLengths) const {
        fastaCollection.populateIdentifierMappingFromFastaCollection(idMap);
        fastaCollection.fillSequenceLengths(sequenceLengths, idMap);
    }
    // Split input fasta files into equal parts and create batches of sequence IDs
    auto getFastaBatches(IdentifierMapping const & idMap) const {
        tsl::hopscotch_map<size_t, std::vector<size_t>> gidToSid;
        for (auto&& elem : idMap.sequenceIDToTuple()) {
            gidToSid[elem.second.gid].emplace_back(elem.first);
        }
        tsl::hopscotch_map<size_t, std::vector<std::vector<size_t>>> gidToBatches;
        for (auto&& elem : gidToSid) {
            for (size_t i = 0; i < config_->batchsize(); ++i) {
                auto chunk = getContainerChunk(elem.second, i, config_->batchsize());
                if (chunk.size()) { gidToBatches[elem.first].emplace_back(chunk); }
            }
        }
        size_t ncombinations = 1;
        std::vector<std::vector<std::vector<size_t>>> batchesPerGid; // cartesianProductByID: T = std::vector<size_t>
        for (auto&& elem : gidToBatches) {
            batchesPerGid.emplace_back(elem.second);
            ncombinations *= elem.second.size();
        }
        std::vector<std::vector<size_t>> batchVector;
        for (size_t i = 0; i < ncombinations; ++i) {
            auto batch = cartesianProductByID<std::vector<size_t>>(i, batchesPerGid);
            std::vector<size_t> flatBatch;
            for (auto&& elem : batch) { flatBatch.insert(flatBatch.end(), elem.begin(), elem.end()); }
            batchVector.emplace_back(flatBatch);
        }
        return batchVector;
    }
    //! Load input fastas from disk into FastaCollection
    auto loadFastas() {
        auto fastaCollection = std::make_shared<FastaCollection>(config_);   // create all FastaRepresentations, including the artificial ones
        std::cout << "[INFO] -- Total number of sequences: " << fastaCollection->numSequences() << std::endl << std::endl;
        output_->outputArtificialSequences(*fastaCollection);
        return fastaCollection;
    }

    bool const allvsall_;
    std::shared_ptr<Configuration> config_;
    std::mutex mutexOutput_;
    std::shared_ptr<Output> output_;
    bool parallelSeedMap_;
    bool const pipelineA_;
};

#endif // SEEDFINDER_H
