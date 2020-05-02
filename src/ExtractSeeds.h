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
#include "GeometricHashing.h"
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
    /*! \param fastaCollection FastaCollection of input fasta files
     * \param seedMap Pointer to a productive child class of SeedMapGeneral, depending on which method to execute
     * \param nThreads Number of threads to use for parallel sequence processing
     * \param noProgressbar Do not print progress bars
     * \param skipMatchCreation Do not create matches after extracting seeds from sequences
     *
     * \details Reads one fasta file at a time, extracting seeds from its sequences in parallel and feeding it into \c seedMap.
     * Calls \c createMatches() from \c seedMap so results are ready after this is done */
    ExtractSeeds(std::shared_ptr<FastaCollection const> fastaCollection,
                 std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType,
                                                TwoBitSeedDataType>> seedMap,
                 int nThreads, bool noProgressbar, bool skipMatchCreation = false);
    virtual ~ExtractSeeds() {}
    //! Run the appropriate extraction
    virtual void extract() { extractAllVsAll(); }
    //! Forwards to \c genome0() method from \c seedMap
    auto const & genome0() { return seedMap_->genome0(); }
    //! Forwards to \c genome1() method from \c seedMap
    auto const & genome1() { return seedMap_->genome1(); }
    //! Forwards to \c numGenomes() method from \c seedMap
    auto numGenomes() const { return seedMap_->numGenomes(); }
    //! Forwards to \c numSequences() method from \c seedMap
    auto numSequences() const { return seedMap_->numSequences(); }

protected:
    //! Factory function for all-vs-all matching
    void extractAllVsAll();
    //! Shared code between \c readSequence() and \c readTuple()
    void extractSeeds(std::string const & sequence,
                      uint32_t sequenceID,
                      std::string const & sequenceName,
                      uint16_t genomeID,
                      std::vector<std::pair<TwoBitKmer<TwoBitKmerDataType>, KmerOccurrence>> & occurrences,
                      size_t span,
                      std::unique_lock<std::mutex> & outputLock);
    //! Process a single sequence, spawned threads execute this function
    void readSequence(std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> seedMapGlobal,
                      typename FastaRepresentation::FastaMapType::const_iterator sequencesIt,
                      typename FastaRepresentation::FastaMapType::const_iterator sequencesEnd,
                      std::string genomeName,
                      ProgressBar & pb);
    //! Helper factory function to check if input data is valid
    void validateInputFiles();

    //! Input sequences
    std::shared_ptr<FastaCollection const> fastaCollection_;
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
    //! Do not create matches after seed extraction if true
    bool skipMatchCreation_;
};



//! Do the complete pipeline (extract seeds, create matches, perform filter) pairwise all-vs-all
/*! This saves memory as the seeds/matches are deleted after each pair was processed.
 * However, this alters the behaviour of the matchLimit filter as this does not work
 * globally anymore. */
template <typename TwoBitKmerDataType, typename TwoBitSeedDataType>
class ExtractSeedsFast : public ExtractSeeds<TwoBitKmerDataType, TwoBitSeedDataType> {
public:
    using SequencePair = std::array<std::reference_wrapper<FastaSequence const>, 2>;
    using SequencePairVector = std::vector<SequencePair>;

    ExtractSeedsFast(std::shared_ptr<Configuration const> config,
                     std::shared_ptr<FastaCollection const> fastaCollection,
                     std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> seedMap,
                     std::shared_ptr<std::ofstream> os)
        : ExtractSeeds<TwoBitKmerDataType, TwoBitSeedDataType>(fastaCollection,
                                                               seedMap,
                                                               config->nThreads(),
                                                               config->quiet(),
                                                               config->performGeometricHashing()),
          config_{config}, mutexOutfile_{}, mutexStdout_{}, os_{os} {}
    //! Run the appropriate extraction
    void extract() { extractPairwise(); }

private:
    void extractPairwise();
    void runPipeline(typename SequencePairVector::const_iterator pairsIt,
                     typename SequencePairVector::const_iterator pairsEnd,
                     ProgressBar & pb);

    //! Shared configuration
    std::shared_ptr<Configuration const> config_;
    //! Lock for json output
    std::mutex mutexOutfile_;
    //! Lock for stdout
    std::mutex mutexStdout_;
    //! output stream for matches
    std::shared_ptr<std::ofstream> os_;
};

#endif // EXTRACTSEEDS_H
