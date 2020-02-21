#ifndef EXTRACTSEEDS_H
#define EXTRACTSEEDS_H

#include <filesystem>
#include <iterator>
#include <memory>
#include <mutex>
#include <thread>

#include "prettyprint/prettyprint.hpp"
#include "Configuration.h"
#include "ContainerChunks.h"
#include "FastaCollection.h"
#include "FastaRepresentation.h"
#include "IdentifierMapping.h"
#include "KmerOccurrence.h"
#include "ProgressBar.h"
#include "SeedMapGeneral.h"
#include "Timestep.h"
#include "TwoBitKmer.h"

//! Reads sequences and extracts seeds from it, feeding seeds into a productive child of SeedMapGeneral and finally calling match creation
template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
class ExtractSeeds {
public:
    //! c'tor (1)
    /*! \param inputFiles Vector of input fasta files
     * \param seedMap Pointer to a productive child class of SeedMapGeneral, depending on which method to execute
     * \param nThreads Number of threads to use for parallel sequence processing
     * \param noProgressbar Do not print progress bars
     *
     * \details Reads one fasta file at a time, extracting seeds from its sequences in parallel and feeding it into \c seedMap.
     * Calls \c createMatches() from \c seedMap so results are ready after this is done */
    ExtractSeeds(std::shared_ptr<FastaCollection const> fastaCollection,
                 std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                TwoBitSeedDataType>> seedMap,
                 int nThreads, bool noProgressbar, bool skipMatchCreation = false);
    //! Forwards to \c genome0() method from \c seedMap
    auto const & genome0() { return seedMap_->genome0(); }
    //! Forwards to \c genome1() method from \c seedMap
    auto const & genome1() { return seedMap_->genome1(); }
    //! Forwards to \c numGenomes() method from \c seedMap
    auto numGenomes() const { return seedMap_->numGenomes(); }
    //! Forwards to \c numSequences() method from \c seedMap
    auto numSequences() const { return seedMap_->numSequences(); }

private:
    //! Factory function for all-vs-all matching
    void extractAllVsAll(std::shared_ptr<FastaCollection const> fastaCollection, bool skipMatchCreation = false);
    //! Shared code between \c readSequence() and \c readTuple()
    void extractSeeds(std::string const & sequence,
                      uint32_t sequenceID,
                      std::string const & sequenceName,
                      uint16_t genomeID,
                      std::string const & genomeName,
                      std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> & seedMapLocal,
                      std::unique_lock<std::mutex> & outputLock);
    //! Process a single sequence, spawned threads execute this function
    void readSequence(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> seedMapGlobal,
                      typename FastaRepresentation::FastaMapType::const_iterator sequencesIt,
                      typename FastaRepresentation::FastaMapType::const_iterator sequencesEnd,
                      std::string genomeName,
                      ProgressBar & pb);
    //! Helper factory function to check if input data is valid
    void validateInputFiles(std::shared_ptr<FastaCollection const> fastaCollection);

    //! Lock for member access of this class
    std::mutex mutexMemberAccess_;
    //! Lock for stdout
    std::mutex mutexOutput_;
    //! Do not print a progress bar if true
    bool noProgressbar_;
    //! Number of cores to use for parallel execution
    int nThreads_;
    //! Feed extracted seeds to this object
    std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                   TwoBitSeedDataType>> seedMap_;
};

#endif // EXTRACTSEEDS_H