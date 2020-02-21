#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <cstddef>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "boost/math/special_functions/binomial.hpp"
#include "boost/program_options.hpp"
#include "json/json.hpp"
#include "prettyprint/prettyprint.hpp"
#include "FastaRepresentation.h"
#include "JsonStream.h"
#include "SpacedSeedMaskCollection.h"
#include "StrongType.h"

namespace fs = std::experimental::filesystem;
namespace po = boost::program_options;

// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
using AllowOverlap = NamedType<bool, struct AllowOverlapTag>;
using ArtificialSequenceSizeFactor = NamedType<size_t, struct ArtificialSequenceSizeFactorTag>;
using CreateAllMatches = NamedType<bool, struct CreateAllMatchesTag>;
using CubeScoreMu = NamedType<double, struct CubeScoreMuTag>;
using CubeScoreThreshold = NamedType<double, struct CubeScoreThresholdTag>;
using DiagonalThreshold = NamedType<double, struct DiagonalThresholdTag>;
using DynamicArtificialSequences = NamedType<bool, struct DynamicArtificialSequencesTag>;
using InputFiles = NamedType<std::vector<std::string>, struct InputFilesTag>;
using Genome1 = NamedType<std::string, struct Genome1Tag>;
using Genome2 = NamedType<std::string, struct Genome2Tag>;
using LocalSearchAreaLength = NamedType<size_t, struct LocalSearchAreaLengthTag>;
using Masks = NamedType<std::vector<std::string>, struct MasksTag>;
using MatchLimit = NamedType<size_t, struct MatchLimitTag>;
using MatchLimitDiscardSeeds = NamedType<bool, struct MatchLimitDiscardSeedsTag>;
using MinMatchDistance = NamedType<size_t, struct MinMatchDistanceTag>;
using NoProgressbar = NamedType<bool, struct NoProgressbarTag>;
using NThreads = NamedType<size_t, struct NThreadsTag>;
using OccurrencePerGenomeMax = NamedType<size_t, struct OccurrencePerGenomeMaxTag>;
using OccurrencePerGenomeMin = NamedType<size_t, struct OccurrencePerGenomeMinTag>;
using OptimalSeed = NamedType<bool, struct OptimalSeedTag>;
using Output = NamedType<fs::path, struct OutputTag>;
using OutputArtificialSequences = NamedType<fs::path, struct OutputArtificialSequencesTag>;
using OutputRunInformation = NamedType<fs::path, struct OutputRunInformationTag>;
using PerformDiagonalFiltering = NamedType<bool, struct PerformDiagonalFilteringTag>;
using PerformGeometricHashing = NamedType<bool, struct PerformGeometricHashingTag>;
using SeedSetSize = NamedType<size_t, struct SeedSetSizeTag>;
using Span = NamedType<size_t, struct SpanTag>;
using TileSize = NamedType<size_t, struct TileSizeTag>;
using Weight = NamedType<size_t, struct WeightTag>;



//! Parses command line parameters and stores the values
class Configuration {
public:
    using ArtificialSequenceLengthsMapType = std::map<std::string, size_t>;
    using MasksType = std::vector<std::string>;

    //! c'tor (1)
    /*! \param argc Argument counter from \c main() function
     * \param argv Argument vector from \c main() function
     *
     * \details Takes care of argument parsing, storage and sanity checking */
    Configuration(int argc, char * argv[]);
    //! c'tor (2)
    /*! \details Explicitly define each member directly */
    Configuration(AllowOverlap allowOverlap,
                  ArtificialSequenceSizeFactor artificialSequenceSizeFactor,
                  CreateAllMatches createAllMatches,
                  CubeScoreMu cubeScoreMu,
                  CubeScoreThreshold cubeScoreThreshold,
                  DiagonalThreshold diagonalThreshold,
                  DynamicArtificialSequences dynamicArtificialSequences,
                  Genome1 genome1,
                  Genome2 genome2,
                  InputFiles inputFiles,
                  LocalSearchAreaLength localSearchAreaLength,
                  Masks masks,
                  MatchLimit matchLimit,
                  MatchLimitDiscardSeeds matchLimitDiscardSeeds,
                  MinMatchDistance minMatchDistance,
                  NoProgressbar noProgressbar,
                  NThreads nThreads,
                  OccurrencePerGenomeMax occurrencePerGenomeMax,
                  OccurrencePerGenomeMin occurrencePerGenomeMin,
                  OptimalSeed optimalSeed,
                  Output output,
                  OutputArtificialSequences outputArtificialSequences,
                  OutputRunInformation outputRunInformation,
                  PerformDiagonalFiltering performDiagonalFiltering,
                  PerformGeometricHashing performGeometricHashing,
                  SeedSetSize seedSetSize,
                  Span span,
                  TileSize tileSize,
                  Weight weight)
        : allowOverlap_{allowOverlap.get()},
          artificialSequenceSizeFactor_{artificialSequenceSizeFactor.get()},
          createAllMatches_{createAllMatches.get()},
          cubeScoreMu_{cubeScoreMu.get()},
          cubeScoreThreshold_{cubeScoreThreshold.get()},
          diagonalThreshold_{diagonalThreshold.get()},
          dynamicArtificialSequences_{dynamicArtificialSequences.get()},
          genome1_{genome1.get()},
          genome2_{genome2.get()},
          inputFiles_{inputFiles.get()},
          localSearchArea_{localSearchAreaLength.get()},
          masks_{masks.get()},
          matchLimit_{matchLimit.get()},
          matchLimitDiscardSeeds_{matchLimitDiscardSeeds.get()},
          minMatchDistance_{minMatchDistance.get()},
          noProgressbar_{noProgressbar.get()},
          nThreads_{nThreads.get()},
          occurrencePerGenomeMax_{occurrencePerGenomeMax.get()},
          occurrencePerGenomeMin_{occurrencePerGenomeMin.get()},
          optimalSeed_{optimalSeed.get()},
          output_{output.get()},
          outputArtificialSequences_{outputArtificialSequences.get()},
          outputRunInformation_{outputRunInformation.get()},
          performDiagonalFiltering_{performDiagonalFiltering.get()},
          performGeometricHashing_{performGeometricHashing.get()},
          seedSetSize_{seedSetSize.get()},
          span_{span.get()},
          tileSize_{tileSize.get()},
          weight_{weight.get()} {}

    //! Getter function for member \c allowOverlap_
    auto allowOverlap() const { return allowOverlap_; }
    //! Getter function for member \c artificialSequenceSizeFactor_
    auto artificialSequenceSizeFactor() const { return artificialSequenceSizeFactor_; }
    //! Flag if only matches from seeds that occur in both genome 0 and 1 should be created
    auto createAllMatches() const { return createAllMatches_; }
    //! Return a JsonValue of a json dict of this configuration with parameters as keys and their respective values
    JsonValue configJson() const;
    //! Getter function for member \c cubeScoreMu_
    auto cubeScoreMu() const { return cubeScoreMu_; }
    //! Getter function for member \c cubeScoreThreshold_
    auto cubeScoreThreshold() const { return cubeScoreThreshold_; }
    //! Getter function for member \c diagonalThreshold_
    auto diagonalThreshold() const { return diagonalThreshold_; }
    //! Getter function for member \c dynamicArtificialSequences_
    auto dynamicArtificialSequences() const { return dynamicArtificialSequences_; }
    //! Getter function for member \c genome1_
    auto const & genome1() const { return genome1_; }
    //! Getter function for member \c genome2_
    auto const & genome2() const { return genome2_; }
    //! Getter function for member \c inputFiles_
    auto const & inputFiles() const { return inputFiles_; }
    //! Getter function for member \c localAreaLength_
    auto localAreaLength() const { return localSearchArea_; }
    //! Getter function for member \c masks_
    auto const & masks() const { return masks_; }
    //! Getter function for member \c matchLimit_
    auto matchLimit() const { return matchLimit_; }
    //! Getter function for member \c matchLimitDiscardSeeds_
    auto matchLimitDiscardSeeds() const { return matchLimitDiscardSeeds_; }
    //! Getter function for member \c minMatchDistance_
    auto minMatchDistance() const { return minMatchDistance_; }
    //! Getter function for member \c nThreads_
    auto nThreads() const { return nThreads_; }
    //! Getter funciton for member \c occurrencePerGenomeMax_
    auto occurrencePerGenomeMax() const { return occurrencePerGenomeMax_; }
    //! Getter funciton for member \c occurrencePerGenomeMin_
    auto occurrencePerGenomeMin() const { return occurrencePerGenomeMin_; }
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
    //! Getter function for member \c noProgressbar_
    auto quiet() const { return noProgressbar_; }
    //! Getter function for member \c seedSetSize_
    auto seedSetSize() const { return seedSetSize_; }
    //! Getter function for member \c span_
    auto span() const { return span_; }
    //! Getter function for member \c tileSize_
    auto tileSize() const { return tileSize_; }
    //! Getter function for member \c weight_
    auto weight() const { return weight_; }
    //! Prints parameter configuration to \c os
    friend std::ostream & operator<<(std::ostream & os, Configuration const & conf);

private:
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
    //! Create artificial sequence(s) of (sum of) lengths of this times the length of the input sequences
    size_t artificialSequenceSizeFactor_;
    //! If true, also create matches from seeds that not occur in genome 0 or 1
    bool createAllMatches_;
    //! [M6] Parameter for Cube scoring
    double cubeScoreMu_;
    //! [M6] Cubes must have at least this score
    double cubeScoreThreshold_;
    //! [M4] Threshold after how many neighbouring matches a match is reported
    double diagonalThreshold_;
    //! Create for each input sequence an artificial sequence of the same length in the respective genome
    bool dynamicArtificialSequences_;
    //! Find matches between this and \c genome2_
    std::string genome1_;
    //! Find matches between this and \c genome1_
    std::string genome2_;
    //! List of input fasta files (one file per genome)
    std::vector<std::string> inputFiles_;
    //! [M4] Length of the area in which to count neighbouring matches
    size_t localSearchArea_;
    //! Optional user defined spaced seed mask(s)
    MasksType masks_;
    //! Max. number of matches (or Link s) to create from a seed
    size_t matchLimit_;
    //! Discard seeds that give more matches (or Link s) than \c matchLimit_ allows
    bool matchLimitDiscardSeeds_;
    //! [M4] Minimal distance between two neighbouring matches
    /*! Has no effect if \c allowOverlap_ is \c true */
    size_t minMatchDistance_;
    //! Don't create progress bars
    bool noProgressbar_;
    //! Number of threads to use in parallel processing steps
    size_t nThreads_;
    //! [M6] At most this many occurrences of a seed in any genome
    size_t occurrencePerGenomeMax_;
    //! [M6] Zero or at least this many occurrences of a seed in any genome
    size_t occurrencePerGenomeMin_;
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
    //! [M3-M4] Number of SpacedSeedMask s to create
    size_t seedSetSize_;
    //! Span (i.e. overall length) of seeds
    size_t span_;
    //! [M6] Tile size for geometricHashing
    size_t tileSize_;
    //! Weight of seeds
    size_t weight_;
};

#endif // CONFIGURATION_H
