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
        allvsall_{false},
        artificialSequenceSizeFactor_(1),
        batchsize_(1),
        createAllMatches_(false),
        cubeLengthCutoff_{30000},
        cubeOutput_{0},
        cubeScoreNormalizationParameter_(30000),
        cubeScoreParameter_(500),
        cubeScoreParameterChunks_(0),
        cubeScoreThreshold_(25),
        diagonalDelta_(5),
        diagonalRho_(200),
        diagonalThreshold_(2),
        dynamicArtificialSequences_(false),
        graphAnnotationFile_(""),
        graphFile_{""},
        genome1_(""),
        genome2_(""),
        hasse_(false),
        inputFiles_(std::vector<std::string>()),
        localSearchAreaLength_(1000),
        maskCollection_{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                         SpacedSeedMaskCollection::Span{5},
                                                                         SpacedSeedMaskCollection::SeedSetSize{1})},
        matchLimit_(10),
        matchLimitDiscardSeeds_(false),
        maxPrefixLength_{12},
        metagraphInterfacePtr_{nullptr},
        minMatchDistance_(0),
        //noProgressbar_(false),
        nThreads_(std::thread::hardware_concurrency()),
        occurrencePerGenomeMax_(ULLONG_MAX),
        occurrencePerGenomeMin_(1),
        occurrencePerSequenceMax_(ULLONG_MAX),
        oldCubeScore_(false),
        optimalSeed_(false),
        output_(""),
        outputArtificialSequences_(""),
        outputRunInformation_(""),
        performDiagonalFiltering_(false),
        performGeometricHashing_(false),
        preAddNeighbouringCubes_{false},
        preHasse_{false},
        preLinkThreshold_{0},
        preMaskCollection_{nullptr},
        preOptimalSeed_{false},
        postSeqential_{false},
        redmask_(false),
        thinning_{5},
        tileSize_(10000),
        verbose_(2),
        yass_{false},
        yassEpsilon_(0.05),
        yassIndel_(0.08),
        yassMutation_(0.15) {
        auto p = (std::thread::hardware_concurrency() > 0)
                    ? std::thread::hardware_concurrency()
                    : 1;
        if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }
        nThreads_ = NThreads(p);
    }

    // create a shared_ptr to a config
    auto makeConfig() const {
        auto config = std::make_shared<Configuration>(allowOverlap_,
                                                      allvsall_,
                                                      artificialSequenceSizeFactor_,
                                                      batchsize_,
                                                      createAllMatches_,
                                                      cubeLengthCutoff_,
                                                      cubeOutput_,
                                                      cubeScoreNormalizationParameter_,
                                                      cubeScoreParameter_,
                                                      cubeScoreParameterChunks_,
                                                      cubeScoreThreshold_,
                                                      diagonalDelta_,
                                                      diagonalRho_,
                                                      diagonalThreshold_,
                                                      dynamicArtificialSequences_,
                                                      genome1_,
                                                      genome2_,
                                                      graphAnnotationFile_,
                                                      graphFile_,
                                                      hasse_,
                                                      inputFiles_,
                                                      localSearchAreaLength_,
                                                      maskCollection_,
                                                      matchLimit_,
                                                      matchLimitDiscardSeeds_,
                                                      maxPrefixLength_,
                                                      metagraphInterfacePtr_,
                                                      minMatchDistance_,
                                                      //noProgressbar_,
                                                      nThreads_,
                                                      occurrencePerGenomeMax_,
                                                      occurrencePerGenomeMin_,
                                                      occurrencePerSequenceMax_,
                                                      oldCubeScore_,
                                                      optimalSeed_,
                                                      output_,
                                                      outputArtificialSequences_,
                                                      outputRunInformation_,
                                                      performDiagonalFiltering_,
                                                      performGeometricHashing_,
                                                      preAddNeighbouringCubes_,
                                                      preHasse_,
                                                      preLinkThreshold_,
                                                      preMaskCollection_,
                                                      preOptimalSeed_,
                                                      postSeqential_,
                                                      redmask_,
                                                      thinning_,
                                                      tileSize_,
                                                      verbose_,
                                                      yass_,
                                                      yassEpsilon_,
                                                      yassIndel_,
                                                      yassMutation_);
        return config;
    }
    void set(AllowOverlap x) {allowOverlap_ = x;}
    void set(Allvsall x) {allvsall_ = x;}
    void set(ArtificialSequenceSizeFactor x) {artificialSequenceSizeFactor_ = x;}
    void set(Batchsize x) {batchsize_ = x;}
    void set(CreateAllMatches x) {createAllMatches_ = x;}
    void set(CubeLengthCutoff x) {cubeLengthCutoff_ = x;}
    void set(CubeOutput x) {cubeOutput_ = x;}
    void set(CubeScoreNormalizationParameter x) {cubeScoreNormalizationParameter_ = x;}
    void set(CubeScoreParameter x) {cubeScoreParameter_ = x;}
    void set(CubeScoreParameterChunks x) {cubeScoreParameterChunks_ = x;}
    void set(CubeScoreThreshold x) {cubeScoreThreshold_ = x;}
    void set(DiagonalDelta x) {diagonalDelta_ = x;}
    void set(DiagonalRho x) {diagonalRho_ = x;}
    void set(DiagonalThreshold x) {diagonalThreshold_ = x;}
    void set(DynamicArtificialSequences x) {dynamicArtificialSequences_ = x;}
    void set(Genome1 x) {genome1_ = x;}
    void set(Genome2 x) {genome2_ = x;}
    void set(GraphAnnotationFile x) {graphAnnotationFile_ = x;}
    void set(GraphFile x) {graphFile_ = x;}
    void set(Hasse x) {hasse_ = x;}
    void set(InputFiles x) {inputFiles_ = x;}
    void set(LocalSearchAreaLength x) {localSearchAreaLength_ = x;}
    void set(MaskCollectionPtr x) {maskCollection_ = x;}
    void set(MatchLimit x) {matchLimit_ = x;}
    void set(MatchLimitDiscardSeeds x) {matchLimitDiscardSeeds_ = x;}
    void set(MaxPrefixLength x) {maxPrefixLength_ = x;}
    void set(MetagraphInterfacePtr x) {metagraphInterfacePtr_ = x;}
    void set(MinMatchDistance x) {minMatchDistance_ = x;}
    //void set(NoProgressbar x) {noProgressbar_ = x;}
    void set(NThreads x) {nThreads_ = x;}
    void set(OccurrencePerGenomeMax x) {occurrencePerGenomeMax_ = x;}
    void set(OccurrencePerGenomeMin x) {occurrencePerGenomeMin_ = x;}
    void set(OccurrencePerSequenceMax x) {occurrencePerSequenceMax_ = x;}
    void set(OldCubeScore x) {oldCubeScore_ = x;}
    void set(OptimalSeed x) {optimalSeed_ = x;}
    void set(OutputPath x) {output_ = x;}
    void set(OutputArtificialSequences x) {outputArtificialSequences_ = x;}
    void set(OutputRunInformation x) {outputRunInformation_ = x;}
    void set(PerformDiagonalFiltering x) {performDiagonalFiltering_ = x;}
    void set(PerformGeometricHashing x) {performGeometricHashing_ = x;}
    void set(PreAddNeighbouringCubes x) {preAddNeighbouringCubes_ = x;}
    void set(PreHasse x) {preHasse_ = x;}
    void set(PreLinkThreshold x) {preLinkThreshold_ = x;}
    void set(PreMaskCollectionPtr x) {preMaskCollection_ = x;}
    void set(PreOptimalSeed x) {preOptimalSeed_ = x;}
    void set(PostSequential x) {postSeqential_ = x;}
    void set(Redmask x) {redmask_ = x;}
    void set(Thinning x) {thinning_ = x;}
    void set(TileSize x) {tileSize_ = x;}
    void set(Verbose x) {verbose_ = x;}
    void set(Yass x) {yass_ = x;}
    void set(YassEpsilon x) {yassEpsilon_ = x;}
    void set(YassIndel x) {yassIndel_ = x;}
    void set(YassMutation x) {yassMutation_ = x;}

    // recursive function to deal with variadic arguments
    template<typename T, typename... Args>
    void set(T t, Args... args) {
        set(t); // set T
        set(args...); // call again, processing the next of args
    }

private:
    AllowOverlap allowOverlap_;
    Allvsall allvsall_;
    ArtificialSequenceSizeFactor artificialSequenceSizeFactor_;
    Batchsize batchsize_;
    CreateAllMatches createAllMatches_;
    CubeLengthCutoff cubeLengthCutoff_;
    CubeOutput cubeOutput_;
    CubeScoreNormalizationParameter cubeScoreNormalizationParameter_;
    CubeScoreParameter cubeScoreParameter_;
    CubeScoreParameterChunks cubeScoreParameterChunks_;
    CubeScoreThreshold cubeScoreThreshold_;
    DiagonalDelta diagonalDelta_;
    DiagonalRho diagonalRho_;
    DiagonalThreshold diagonalThreshold_;
    DynamicArtificialSequences dynamicArtificialSequences_;
    GraphAnnotationFile graphAnnotationFile_;
    GraphFile graphFile_;
    Genome1 genome1_;
    Genome2 genome2_;
    Hasse hasse_;
    InputFiles inputFiles_;
    LocalSearchAreaLength localSearchAreaLength_;
    MaskCollectionPtr maskCollection_;
    MatchLimit matchLimit_;
    MatchLimitDiscardSeeds matchLimitDiscardSeeds_;
    MaxPrefixLength maxPrefixLength_;
    MetagraphInterfacePtr metagraphInterfacePtr_;
    MinMatchDistance minMatchDistance_;
    //NoProgressbar noProgressbar_;
    NThreads nThreads_;
    OccurrencePerGenomeMax occurrencePerGenomeMax_;
    OccurrencePerGenomeMin occurrencePerGenomeMin_;
    OccurrencePerSequenceMax occurrencePerSequenceMax_;
    OldCubeScore oldCubeScore_;
    OptimalSeed optimalSeed_;
    OutputPath output_;
    OutputArtificialSequences outputArtificialSequences_;
    OutputRunInformation outputRunInformation_;
    PerformDiagonalFiltering performDiagonalFiltering_;
    PerformGeometricHashing performGeometricHashing_;
    PreAddNeighbouringCubes preAddNeighbouringCubes_;
    PreHasse preHasse_;
    PreLinkThreshold preLinkThreshold_;
    PreMaskCollectionPtr preMaskCollection_;
    PreOptimalSeed preOptimalSeed_;
    PostSequential postSeqential_;
    Redmask redmask_;
    Thinning thinning_;
    TileSize tileSize_;
    Verbose verbose_;
    Yass yass_;
    YassEpsilon yassEpsilon_;
    YassIndel yassIndel_;
    YassMutation yassMutation_;
};



// variadic function to get a custom config with only the neccessary adjustments
template<typename T, typename... Args>
auto customConfiguration(T t, Args... args) {
    auto configBuilder = ConfigBuilder();
    configBuilder.set(t, args...);
    return configBuilder.makeConfig();
}

#endif // CONFIGURATIONGENERATOR_H
