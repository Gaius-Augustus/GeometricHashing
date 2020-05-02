#ifndef CONFIGURATIONGENERATOR_H
#define CONFIGURATIONGENERATOR_H

#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "Configuration.h"

// https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
// split strings at each space and insert parts into argvStr
inline void addParameterFromString(std::vector<std::string> & argvStr, std::string const & str) {
    std::istringstream iss(str);
    argvStr.insert(argvStr.end(),
                   std::istream_iterator<std::string>(iss),
                   std::istream_iterator<std::string>());
}



inline auto configurationPtrFromParameters(std::vector<std::string> const & argvStr) {
    // convert into C++-style argv
    char** argv = new char*[argvStr.size() + 1];
    for (size_t i = 0; i < argvStr.size(); ++i) {
        argv[i] = new char[argvStr.at(i).size() + 1];   // +1 for null termination
        argvStr.at(i).copy(argv[i], argvStr.at(i).size());
        argv[i][argvStr.at(i).size()] = 0;  // null terminating
    }
    argv[argvStr.size()] = nullptr;

    //for (size_t i = 0; i < argvStr.size(); ++i) { std::cout << argv[i] << std::endl; }

    auto config = std::make_shared<Configuration>(argvStr.size(), argv);

    // cleanup raw pointer mess
    for (size_t i = 0; i < argvStr.size() + 1; ++i) { delete[] argv[i]; }
    delete[] argv;

    return config;
}



inline auto configurationFromParameters(std::vector<std::string> const & argvStr) {
    auto config = configurationPtrFromParameters(argvStr);
    return *config;
}



inline auto generateConfiguration(std::string const & str) {
    std::vector<std::string> argvStr;
    addParameterFromString(argvStr, str);
    return configurationPtrFromParameters(argvStr);
}

/* Builder class that stores default values for each configuration item
     Each item can be modified before building a config */
class ConfigBuilder {
public:
    // c'tor -- set defaults, these are used if not changed by set()
    ConfigBuilder() :
        allowOverlap_(false),
        artificialSequenceSizeFactor_(1),
        createAllMatches_(false),
        cubeAreaCutoff_{300000000},
        cubeScoreNormalizationParameter_(300000000),
        cubeScoreParameter_(500),
        cubeScoreThreshold_(25),
        diagonalThreshold_(2),
        dynamicArtificialSequences_(false),
        fast_(false),
        fastBatchsize_(0),
        genome1_(""),
        genome2_(""),
        inputFiles_(std::vector<std::string>()),
        localSearchAreaLength_(1000),
        masks_(std::vector<std::string>()),
        matchLimit_(10),
        matchLimitDiscardSeeds_(false),
        minMatchDistance_(0),
        noProgressbar_(false),
        nThreads_(std::thread::hardware_concurrency()),
        occurrencePerGenomeMax_(ULLONG_MAX),
        occurrencePerGenomeMin_(1),
        oldCubeScore_(false),
        optimalSeed_(false),
        output_(""),
        outputArtificialSequences_(""),
        outputRunInformation_(""),
        performDiagonalFiltering_(false),
        performGeometricHashing_(false),
        seedSetSize_(1),
        span_(5),
        tileSize_(10000),
        weight_(5) {
        auto p = (std::thread::hardware_concurrency() > 0)
                    ? std::thread::hardware_concurrency()
                    : 1;
        if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }
        nThreads_ = NThreads(p);
    }

    // create a shared_ptr to a config
    auto makeConfig() const {
        auto config = std::make_shared<Configuration>(allowOverlap_,
                                                      artificialSequenceSizeFactor_,
                                                      createAllMatches_,
                                                      cubeAreaCutoff_,
                                                      cubeScoreNormalizationParameter_,
                                                      cubeScoreParameter_,
                                                      cubeScoreThreshold_,
                                                      diagonalThreshold_,
                                                      dynamicArtificialSequences_,
                                                      fast_,
                                                      fastBatchsize_,
                                                      genome1_,
                                                      genome2_,
                                                      inputFiles_,
                                                      localSearchAreaLength_,
                                                      masks_,
                                                      matchLimit_,
                                                      matchLimitDiscardSeeds_,
                                                      minMatchDistance_,
                                                      noProgressbar_,
                                                      nThreads_,
                                                      occurrencePerGenomeMax_,
                                                      occurrencePerGenomeMin_,
                                                      oldCubeScore_,
                                                      optimalSeed_,
                                                      output_,
                                                      outputArtificialSequences_,
                                                      outputRunInformation_,
                                                      performDiagonalFiltering_,
                                                      performGeometricHashing_,
                                                      seedSetSize_,
                                                      span_,
                                                      tileSize_,
                                                      weight_);
        return config;
    }
    void set(AllowOverlap x) {allowOverlap_ = x;}
    void set(ArtificialSequenceSizeFactor x) {artificialSequenceSizeFactor_ = x;}
    void set(CreateAllMatches x) {createAllMatches_ = x;}
    void set(CubeAreaCutoff x) {cubeAreaCutoff_ = x;}
    void set(CubeScoreNormalizationParameter x) {cubeScoreNormalizationParameter_ = x;}
    void set(CubeScoreParameter x) {cubeScoreParameter_ = x;}
    void set(CubeScoreThreshold x) {cubeScoreThreshold_ = x;}
    void set(DiagonalThreshold x) {diagonalThreshold_ = x;}
    void set(DynamicArtificialSequences x) {dynamicArtificialSequences_ = x;}
    void set(Fast x) {fast_ = x;}
    void set(FastBatchsize x) {fastBatchsize_ = x;}
    void set(Genome1 x) {genome1_ = x;}
    void set(Genome2 x) {genome2_ = x;}
    void set(InputFiles x) {inputFiles_ = x;}
    void set(LocalSearchAreaLength x) {localSearchAreaLength_ = x;}
    void set(Masks x) {masks_ = x;}
    void set(MatchLimit x) {matchLimit_ = x;}
    void set(MatchLimitDiscardSeeds x) {matchLimitDiscardSeeds_ = x;}
    void set(MinMatchDistance x) {minMatchDistance_ = x;}
    void set(NoProgressbar x) {noProgressbar_ = x;}
    void set(NThreads x) {nThreads_ = x;}
    void set(OccurrencePerGenomeMax x) {occurrencePerGenomeMax_ = x;}
    void set(OccurrencePerGenomeMin x) {occurrencePerGenomeMin_ = x;}
    void set(OldCubeScore x) {oldCubeScore_ = x;}
    void set(OptimalSeed x) {optimalSeed_ = x;}
    void set(Output x) {output_ = x;}
    void set(OutputArtificialSequences x) {outputArtificialSequences_ = x;}
    void set(OutputRunInformation x) {outputRunInformation_ = x;}
    void set(PerformDiagonalFiltering x) {performDiagonalFiltering_ = x;}
    void set(PerformGeometricHashing x) {performGeometricHashing_ = x;}
    void set(SeedSetSize x) {seedSetSize_ = x;}
    void set(Span x) {span_ = x;}
    void set(TileSize x) {tileSize_ = x;}
    void set(Weight x) {weight_ = x;}
    // recursive function to deal with variadic arguments
    template<typename T, typename... Args>
    void set(T t, Args... args) {
        set(t); // set T
        set(args...); // call again, processing the next of args
    }

private:
    AllowOverlap allowOverlap_;
    ArtificialSequenceSizeFactor artificialSequenceSizeFactor_;
    CreateAllMatches createAllMatches_;
    CubeAreaCutoff cubeAreaCutoff_;
    CubeScoreNormalizationParameter cubeScoreNormalizationParameter_;
    CubeScoreParameter cubeScoreParameter_;
    CubeScoreThreshold cubeScoreThreshold_;
    DiagonalThreshold diagonalThreshold_;
    DynamicArtificialSequences dynamicArtificialSequences_;
    Fast fast_;
    FastBatchsize fastBatchsize_;
    Genome1 genome1_;
    Genome2 genome2_;
    InputFiles inputFiles_;
    LocalSearchAreaLength localSearchAreaLength_;
    Masks masks_;
    MatchLimit matchLimit_;
    MatchLimitDiscardSeeds matchLimitDiscardSeeds_;
    MinMatchDistance minMatchDistance_;
    NoProgressbar noProgressbar_;
    NThreads nThreads_;
    OccurrencePerGenomeMax occurrencePerGenomeMax_;
    OccurrencePerGenomeMin occurrencePerGenomeMin_;
    OldCubeScore oldCubeScore_;
    OptimalSeed optimalSeed_;
    Output output_;
    OutputArtificialSequences outputArtificialSequences_;
    OutputRunInformation outputRunInformation_;
    PerformDiagonalFiltering performDiagonalFiltering_;
    PerformGeometricHashing performGeometricHashing_;
    SeedSetSize seedSetSize_;
    Span span_;
    TileSize tileSize_;
    Weight weight_;
};



// variadic function to get a custom config with only the neccessary adjustments
template<typename T, typename... Args>
auto customConfiguration(T t, Args... args) {
    auto configBuilder = ConfigBuilder();
    configBuilder.set(t, args...);
    return configBuilder.makeConfig();
}

#endif // CONFIGURATIONGENERATOR_H
