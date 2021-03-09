#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "boost/math/special_functions/binomial.hpp"
#include "boost/program_options.hpp"
#include "mabl3/Timestep.h"
#include "computeYassParameters.h"
#include "MetagraphInterface.h"
#include "nlohmann/json.hpp"
#include "prettyprint.hpp"
#include "FastaRepresentation.h"
#include "mabl3/JsonStream.h"
#include "SpacedSeedMaskCollection.h"
#include "StrongType.h"

namespace fs = std::filesystem;
namespace po = boost::program_options;
using namespace mabl3;



// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
using AllowOverlap = NamedType<bool, struct AllowOverlapTag>;
using Allvsall = NamedType<bool, struct AllvsallTag>;
using ArtificialSequenceSizeFactor = NamedType<size_t, struct ArtificialSequenceSizeFactorTag>;
using Batchsize = NamedType<size_t, struct BatchsizeTag>;
using CreateAllMatches = NamedType<bool, struct CreateAllMatchesTag>;
using CubeLengthCutoff = NamedType<size_t, struct CubeLengthCutoffTag>;
using CubeOutput = NamedType<size_t, struct CubeOutputTag>;
using CubeScoreNormalizationParameter = NamedType<size_t, struct CubeScoreNormalizationParameterTag>;
using CubeScoreParameter = NamedType<size_t, struct CubeScoreParameterTag>;
using CubeScoreParameterChunks = NamedType<size_t, struct CubeScoreParameterChunksTag>;
using CubeScoreThreshold = NamedType<double, struct CubeScoreThresholdTag>;
using DiagonalDelta= NamedType<size_t, struct DiagonalDeltaTag>;
using DiagonalRho = NamedType<size_t, struct DiagonalRhoTag>;
using DiagonalThreshold = NamedType<double, struct DiagonalThresholdTag>;
using DynamicArtificialSequences = NamedType<bool, struct DynamicArtificialSequencesTag>;
using InputFiles = NamedType<std::vector<std::string>, struct InputFilesTag>;
using Genome1 = NamedType<std::string, struct Genome1Tag>;
using Genome2 = NamedType<std::string, struct Genome2Tag>;
using GraphAnnotationFile = NamedType<fs::path, struct GraphAnnotationFileTag>;
using GraphFile = NamedType<fs::path, struct GraphFileTag>;
using Hasse = NamedType<bool, struct HasseTag>;
using LocalSearchAreaLength = NamedType<size_t, struct LocalSearchAreaLengthTag>;
using MaskCollectionPtr = NamedType<std::shared_ptr<SpacedSeedMaskCollection const>, struct MaskCollectionPtrTag>;
using MatchLimit = NamedType<size_t, struct MatchLimitTag>;
using MatchLimitDiscardSeeds = NamedType<bool, struct MatchLimitDiscardSeedsTag>;
using MaxPrefixLength = NamedType<size_t, struct MaxPrefixLengthTag>;
using MetagraphInterfacePtr = NamedType<std::shared_ptr<MetagraphInterface const>, struct MetagraphInterfacePtrTag>;
using MinMatchDistance = NamedType<size_t, struct MinMatchDistanceTag>;
using NThreads = NamedType<size_t, struct NThreadsTag>;
using OccurrencePerGenomeMax = NamedType<size_t, struct OccurrencePerGenomeMaxTag>;
using OccurrencePerGenomeMin = NamedType<size_t, struct OccurrencePerGenomeMinTag>;
using OccurrencePerSequenceMax = NamedType<size_t, struct OccurrencePerSequenceMaxTag>;
using OldCubeScore = NamedType<bool, struct OldCubeScoreTag>;
using OptimalSeed = NamedType<bool, struct OptimalSeedTag>;
using OutputPath = NamedType<fs::path, struct OutputTag>;
using OutputArtificialSequences = NamedType<fs::path, struct OutputArtificialSequencesTag>;
using OutputRunInformation = NamedType<fs::path, struct OutputRunInformationTag>;
using PerformDiagonalFiltering = NamedType<bool, struct PerformDiagonalFilteringTag>;
using PerformGeometricHashing = NamedType<bool, struct PerformGeometricHashingTag>;
using PreAddNeighbouringCubes = NamedType<bool, struct PreAddNeighbouringCubesTag>;
using PreHasse = NamedType<bool, struct PreHasseTag>;
using PostSequential = NamedType<bool, struct PostSeqentialTag>;
using PreLinkThreshold = NamedType<size_t, struct PreLinkThresholdTag>;
using PreMaskCollectionPtr = NamedType<std::shared_ptr<SpacedSeedMaskCollection const>, struct PreMaskCollectionPtrTag>;
using PreOptimalSeed = NamedType<bool, struct PreOptimalSeedTag>;
using Redmask = NamedType<bool, struct RedmaksTag>;
using Thinning = NamedType<size_t, struct ThinningTag>;
using TileSize = NamedType<size_t, struct TileSizeTag>;
using Verbose = NamedType<size_t, struct VerboseTag>;
using Yass = NamedType<bool, struct YassTag>;
using YassEpsilon = NamedType<double, struct YassEpsilonTag>;
using YassIndel = NamedType<double, struct YassIndelTag>;
using YassMutation = NamedType<double, struct YassMutationTag>;



//! Parses command line parameters and stores the values
class Configuration {
public:
    using ArtificialSequenceLengthsMapType = std::map<std::string, size_t>;
    //using MasksType = std::vector<std::string>;

    //! c'tor (1)
    /*! \param argc Argument counter from \c main() function
     * \param argv Argument vector from \c main() function
     *
     * \details Takes care of argument parsing, storage and sanity checking */
    Configuration(int argc, char * argv[]);
    //! c'tor (2)
    /*! \details Explicitly define each member directly */
    Configuration(AllowOverlap allowOverlap,
                  Allvsall allvsall,
                  ArtificialSequenceSizeFactor artificialSequenceSizeFactor,
                  Batchsize batchsize,
                  CreateAllMatches createAllMatches,
                  CubeLengthCutoff cubeLengthCutoff,
                  CubeOutput cubeOutput,
                  CubeScoreNormalizationParameter cubeScoreNormalizationParameter,
                  CubeScoreParameter cubeScoreParameter,
                  CubeScoreParameterChunks cubeScoreParameterChunks,
                  CubeScoreThreshold cubeScoreThreshold,
                  DiagonalDelta diagonalDelta,
                  DiagonalRho diagonalRho,
                  DiagonalThreshold diagonalThreshold,
                  DynamicArtificialSequences dynamicArtificialSequences,
                  Genome1 genome1,
                  Genome2 genome2,
                  GraphAnnotationFile graphAnnotationFile,
                  GraphFile graphFile,
                  Hasse hasse,
                  InputFiles inputFiles,
                  LocalSearchAreaLength localSearchAreaLength,
                  MaskCollectionPtr maskCollection,
                  MatchLimit matchLimit,
                  MatchLimitDiscardSeeds matchLimitDiscardSeeds,
                  MaxPrefixLength maxPrefixLength,
                  MetagraphInterfacePtr metagraphInterface,
                  MinMatchDistance minMatchDistance,
                  NThreads nThreads,
                  OccurrencePerGenomeMax occurrencePerGenomeMax,
                  OccurrencePerGenomeMin occurrencePerGenomeMin,
                  OccurrencePerSequenceMax occurrencePerSequenceMax,
                  OldCubeScore oldCubeScore,
                  OptimalSeed optimalSeed,
                  OutputPath output,
                  OutputArtificialSequences outputArtificialSequences,
                  OutputRunInformation outputRunInformation,
                  PerformDiagonalFiltering performDiagonalFiltering,
                  PerformGeometricHashing performGeometricHashing,
                  PreAddNeighbouringCubes preAddNeighbouringCubes,
                  PreHasse preHasse,
                  PreLinkThreshold preLinkThreshold,
                  PreMaskCollectionPtr preMaskCollection,
                  PreOptimalSeed preOptimalSeed,
                  PostSequential postSequential,
                  Redmask redmask,
                  Thinning thinning,
                  TileSize tileSize,
                  Verbose verbose,
                  Yass yass,
                  YassEpsilon yassEpsilon,
                  YassIndel yassIndel,
                  YassMutation yassMutation)
        : allowOverlap_{allowOverlap.get()},
          allvsall_{allvsall.get()},
          artificialSequenceSizeFactor_{artificialSequenceSizeFactor.get()},
          batchsize_{batchsize.get()},
          createAllMatches_{createAllMatches.get()},
          cubeLengthCutoff_{cubeLengthCutoff.get()},
          cubeOutput_{cubeOutput.get()},
          cubeScoreNormalizationParameter_{cubeScoreNormalizationParameter.get()},
          cubeScoreParameter_{cubeScoreParameter.get()},
          cubeScoreParameterChunks_{cubeScoreParameterChunks.get()},
          cubeScoreThreshold_{cubeScoreThreshold.get()},
          diagonalDelta_{diagonalDelta.get()},
          diagonalRho_{diagonalRho.get()},
          diagonalThreshold_{diagonalThreshold.get()},
          dynamicArtificialSequences_{dynamicArtificialSequences.get()},
          genome1_{genome1.get()},
          genome2_{genome2.get()},
          graphAnnotationFile_{graphAnnotationFile.get()},
          graphFile_{graphFile.get()},
          hasse_{hasse.get()},
          inputFiles_{inputFiles.get()},
          localSearchArea_{localSearchAreaLength.get()},
          maskCollection_{maskCollection.get()},
          matchLimit_{matchLimit.get()},
          matchLimitDiscardSeeds_{matchLimitDiscardSeeds.get()},
          maxPrefixLength_{maxPrefixLength.get()},
          metagraphInterface_{metagraphInterface.get()},
          minMatchDistance_{minMatchDistance.get()},
          nThreads_{nThreads.get()},
          occurrencePerGenomeMax_{occurrencePerGenomeMax.get()},
          occurrencePerGenomeMin_{occurrencePerGenomeMin.get()},
          occurrencePerSequenceMax_{occurrencePerSequenceMax.get()},
          oldCubeScore_{oldCubeScore.get()},
          optimalSeed_{optimalSeed.get()},
          output_{output.get()},
          outputArtificialSequences_{outputArtificialSequences.get()},
          outputRunInformation_{outputRunInformation.get()},
          performDiagonalFiltering_{performDiagonalFiltering.get()},
          performGeometricHashing_{performGeometricHashing.get()},
          preAddNeighbouringCubes_{preAddNeighbouringCubes.get()},
          preHasse_{preHasse.get()},
          preLinkThreshold_{preLinkThreshold.get()},
          preMaskCollection_{preMaskCollection.get()},
          preOptimalSeed_{preOptimalSeed.get()},
          postSequential_{postSequential.get()},
          redmask_{redmask.get()},
          thinning_{thinning.get()},
          tileSize_{tileSize.get()},
          verbose_{verbose.get()},
          yass_{yass.get()},
          yassEpsilon_{yassEpsilon.get()},
          yassIndel_{yassIndel.get()},
          yassMutation_{yassMutation.get()} {}

    //! Getter function for member \c allowOverlap_
    auto allowOverlap() const { return allowOverlap_; }
    //! Getter function for member \c allvsall_
    auto allvsall() const { return allvsall_; }
    //! Getter function for member \c artificialSequenceSizeFactor_
    auto artificialSequenceSizeFactor() const { return artificialSequenceSizeFactor_; }
    //! Getter function for member \c batchsize_
    auto const & batchsize() const { return batchsize_; }
    //! Flag if only matches from seeds that occur in both genome 0 and 1 should be created
    auto createAllMatches() const { return createAllMatches_; }
    //! Return a JsonValue of a json dict of this configuration with parameters as keys and their respective values
    JsonValue configJson() const;
    //! Getter function for member \c cubeAreaCutoff_
    auto cubeLengthCutoff() const { return cubeLengthCutoff_; }
    //! Getter function for member \c cubeOutput_
    auto cubeOutput() const  { return cubeOutput_; }
    //! Getter function for member \c cubeScoreNormalizationParameter_
    auto cubeScoreNormalizationParameter() const { return cubeScoreNormalizationParameter_; }
    //! Getter function for member \c cubeScoreParameter_
    auto cubeScoreParameter() const { return cubeScoreParameter_; }
    //! Getter function for member \c cubeScoreParameterChunks_
    auto cubeScoreParameterChunks() const { return cubeScoreParameterChunks_; }
    //! Getter function for member \c cubeScoreThreshold_
    auto cubeScoreThreshold() const { return cubeScoreThreshold_; }
    //! Getter function for member \c diagonalThreshold_
    auto diagonalThreshold() const { return diagonalThreshold_; }
    //! Getter function for member \c deltaDiagonal_
    auto diagonalDelta() const { return diagonalDelta_; }
    //! Getter function for member \c deltaRho_
    auto diagonalRho() const { return diagonalRho_; }
    //! Getter function for member \c dynamicArtificialSequences_
    auto dynamicArtificialSequences() const { return dynamicArtificialSequences_; }
    //! Getter function for member \c genome1_
    auto const & genome1() const { return genome1_; }
    //! Getter function for member \c genome2_
    auto const & genome2() const { return genome2_; }
    //! Getter function for member \c graphAnnotationFile_
    auto const & graphAnnotationFile() const { return graphAnnotationFile_; }
    //! Getter function for member \c graphFile_
    auto const & graphFile() const { return graphFile_; }
    //! Getter function for member \c hasse_
    auto hasse() const { return hasse_; }
    //! Getter function for member \c inputFiles_
    auto const & inputFiles() const { return inputFiles_; }
    //! Getter function for member \c localAreaLength_
    auto localAreaLength() const { return localSearchArea_; }
    //! Getter function for member \c masks_
    auto masks() const { return maskCollection_->masksAsString(); }
    //! Getter function for member \c maskCollection_
    auto maskCollection() const { return maskCollection_; }
    //! Getter function for member \c matchLimit_
    auto matchLimit() const { return matchLimit_; }
    //! Getter function for member \c matchLimitDiscardSeeds_
    auto matchLimitDiscardSeeds() const { return matchLimitDiscardSeeds_; }
    //! Getter function for member \c maxPrefixLength_
    auto maxPrefixLength() const { return maxPrefixLength_; }
    //! Getter function for member \c metagraphInterface_
    auto metagraphInterface() const { return metagraphInterface_; }
    //! Getter function for member \c minMatchDistance_
    auto minMatchDistance() const { return minMatchDistance_; }
    //! Getter function for member \c nThreads_
    auto nThreads() const { return nThreads_; }
    //! Getter funciton for member \c occurrencePerGenomeMax_
    auto occurrencePerGenomeMax() const { return occurrencePerGenomeMax_; }
    //! Getter funciton for member \c occurrencePerGenomeMin_
    auto occurrencePerGenomeMin() const { return occurrencePerGenomeMin_; }
    //! Getter function for member \c occurrencePerSequenceMax_
    auto occurrencePerSequenceMax() const { return occurrencePerSequenceMax_; }
    //! Getter function for member \c oldCubeScore_
    auto oldCubeScore() const { return oldCubeScore_; }
    //! Getter function for member \c optimalSeed_
    auto optimalSeed() const { return optimalSeed_; }
    //! Getter function for member \c output_
    auto const & output() const { return output_; }
    //! Getter function for member \c outputArtificialSequences_
    auto const & outputArtificialSequences() const { return outputArtificialSequences_; }
    //! Getter function for member \c outputRunInformation_
    auto const & outputRunInformation() const { return outputRunInformation_; }
    //! Getter fucntion for member \c performDiagonalFiltering_
    auto performDiagonalFiltering() const { return performDiagonalFiltering_; }
    //! Getter fucntion for member \c performGeometricHashing_
    auto performGeometricHashing() const { return performGeometricHashing_; }
    //! Getter function for member \c preAddNeighbouringCubes_
    auto preAddNeighbouringCubes() const { return preAddNeighbouringCubes_; }
    //! Getter function for member \c preHasse_
    auto preHasse() const { return preHasse_; }
    //! Getter function for member \c preLinkThreshold_
    auto preLinkThreshold() const { return preLinkThreshold_; }
    //! Getter function for member \c preMaskCollection_
    auto preMaskCollection() const { return preMaskCollection_; }
    //! Getter function for member \c preOptimalSeed_
    auto preOptimalSeed() const { return preOptimalSeed_; }
    //! Getter function for member \c postSequential_
    auto postSequential() const { return postSequential_; }
    //! Getter function for member \c redmask_
    auto redmask() const { return redmask_; }
    //! Forward to getter function for size of SpacedSeedMaskCollection
    auto seedSetSize() const { return maskCollection_->size(); }
    //! Forward to getter function for maxSpan of SpacedSeedMaskCollection
    auto span() const { return maskCollection_->maxSpan(); }
    //! Getter function for member \c thinning_
    auto thinning() const { return thinning_; }
    //! Getter function for member \c tileSize_
    auto tileSize() const { return tileSize_; }
    //! Getter function for member \c verbose_
    auto verbose() const { return verbose_; }
    //! Forward to getter function for weight  of SpacedSeedMaskCollection
    auto weight() const { return maskCollection_->weight(); }
    //! Getter function for member \c yass_
    auto yass() const { return yass_; }
    //! Getter function for member \c yassEpsilon_
    auto yassEpsilon() const { return yassEpsilon_; }
    //! Getter function for member \c yassIndel_
    auto yassIndel() const { return yassIndel_; }
    //! Getter function for member \c yassMutation_
    auto yassMutation() const { return yassMutation_; }
    //! Prints parameter configuration to \c os
    friend std::ostream & operator<<(std::ostream & os, Configuration const & conf);

protected:
    //! Check if all CLI options are valid and update member variables accordingly
    void validateOptions(po::variables_map & vm);
    //! Check if a user option is in a certain value range and possibly cast to desired datatype
    /*! Make sure that implicit and explicit cast between \c InputType and \c OutputType is possible */
    template <typename InputType, typename OutputType>
    auto castWithBoundaryCheck(po::variables_map & vm, std::string const & key, InputType min, InputType max) {
        auto i = vm[key].as<InputType>();
        if (i < min || max < i) { throw std::runtime_error("[ERROR] -- '--" + key + "' must be in [" + std::to_string(min) + ", " + std::to_string(max) +"]"); }
        return static_cast<OutputType>(i);
    }


    //! [M4] Allow neighbouring matches to overlap (does not add overlapping fractions to match count)
    /*! \c minMatchDistance_ has no effect if this is \c true */
    bool allowOverlap_;
    //! Create Links from all-vs-all sequences (as opposed to one refseq against all non-ref at a time)
    bool allvsall_;
    //! Create artificial sequence(s) of (sum of) lengths of this times the length of the input sequences
    size_t artificialSequenceSizeFactor_;
    //! Run pipeline from batches of this size
    size_t batchsize_;
    //! If true, also create matches from seeds that not occur in genome 0 or 1
    bool createAllMatches_;
    //! [M6] Parameter for cube score computation
    size_t cubeLengthCutoff_;
    //! [M6] Output cubes to file
    size_t cubeOutput_;
    //! [M6] Parameter for Cube score normalization
    size_t cubeScoreNormalizationParameter_;
    //! [M6] Parameter for Cube scoring
    size_t cubeScoreParameter_;
    //! [M6] Parameter 2 for Cube scoring
    size_t cubeScoreParameterChunks_;
    //! [M6] Cubes must have at least this score
    double cubeScoreThreshold_;
    //! [M4] Diagonal difference parameter
    size_t diagonalDelta_;
    //! [M4] Position difference parameter
    size_t diagonalRho_;
    //! [M4] Threshold after how many neighbouring matches a match is reported
    double diagonalThreshold_;
    //! Create for each input sequence an artificial sequence of the same length in the respective genome
    bool dynamicArtificialSequences_;
    //! Find matches between this and \c genome2_
    std::string genome1_;
    //! Find matches between this and \c genome1_
    std::string genome2_;
    //! Path to metagraph annotation
    fs::path graphAnnotationFile_;
    //! Path to metagraph
    fs::path graphFile_;
    //! [M6] perform hasse subcube stuff
    bool hasse_;
    //! List of input fasta files (one file per genome)
    std::vector<std::string> inputFiles_;
    //! [M4] Length of the area in which to count neighbouring matches
    size_t localSearchArea_;
    //! Mask collection for use in seed finding
    std::shared_ptr<SpacedSeedMaskCollection const> maskCollection_;
    //! Max. number of matches (or Link s) to create from a seed
    size_t matchLimit_;
    //! Discard seeds that give more matches (or Link s) than \c matchLimit_ allows
    bool matchLimitDiscardSeeds_;
    //! For SeedOccurrenceMap, determines base RAM usage vs. query time
    size_t maxPrefixLength_;
    //! Interface to a metagraph instance
    std::shared_ptr<MetagraphInterface const> metagraphInterface_;
    //! [M4] Minimal distance between two neighbouring matches
    /*! Has no effect if \c allowOverlap_ is \c true */
    size_t minMatchDistance_;
    //! Number of threads to use in parallel processing steps
    size_t nThreads_;
    //! At most this many occurrences of a seed in any genome
    size_t occurrencePerGenomeMax_;
    //! Zero or at least this many occurrences of a seed in any genome
    size_t occurrencePerGenomeMin_;
    //! At most this many occurrences of a seed in a sequence for the seq. to be considered in link creation
    size_t occurrencePerSequenceMax_;
    //! Use old scoring algorithm for GH
    bool oldCubeScore_;
    //! Use pre-computed optimal seed instead of random one
    bool optimalSeed_;
    //! File to store json results in
    fs::path output_;
    //! File to store generated artificial sequences
    fs::path outputArtificialSequences_;
    //! File to store run info (json)
    fs::path outputRunInformation_;
    //! Flag to perform M4
    bool performDiagonalFiltering_;
    //! Flag to perform M6
    bool performGeometricHashing_;
    //! If set, for all found relevant cubes all the neighbours are added as well in pre-GH
    bool preAddNeighbouringCubes_;
    //! [M5] If set, allow incomplete cubes in pre-filtering step (when >= 3 input genomes)
    bool preHasse_;
    //! [M5] Link count threshold for pre-filtering step
    size_t preLinkThreshold_;
    //! [M5] If set, perform pre-filtering step in GH
    std::shared_ptr<SpacedSeedMaskCollection const> preMaskCollection_;
    //! [M5] Optimal seeds used in pre-filtering step
    bool preOptimalSeed_;
    //! [M5] If set, run the second GH step sequentially rather than in parallel
    bool postSequential_;
    //! Discard low-complexity seeds (only one or two nt in seed), like YASS
    bool redmask_;
    //! Discard roughly 1/thinning_ of input k-mers
    size_t thinning_;
    //! [M6] Tile size for geometricHashing
    size_t tileSize_;
    //! Detail level of output (0 - none, 1 - static, 2 - static + progress bars
    size_t verbose_;
    //! [M4] perform yass-like
    bool yass_;
    //! [M4] epsilon parameter for yass-like filtering
    double yassEpsilon_;
    //! [M4] indel rate parameter for yass-like filtering
    double yassIndel_;
    //! [M4] mutation rate parameter for yass-like filtering
    double yassMutation_;
};

#endif // CONFIGURATION_H
