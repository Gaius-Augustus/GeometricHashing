#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "catch2/catch.hpp"
#include "../Configuration.h"
#include "ConfigurationGenerator.h"
#include "WriteTestdata.h"



auto throwsOnConstruction(std::vector<std::string> const & argvStr) {
    auto thrown = false;
    try {
        configurationFromParameters(argvStr);
    } catch (std::exception const & e) {
        thrown = true;
        std::cout << "[INFO] -- throwsOnConstruction -- Message: " << e.what() << std::endl;
    }
    return thrown;
}



TEST_CASE("Configuration") {
    std::vector<std::string> argvStr;
    argvStr.emplace_back("seedFinding");    // program call
    argvStr.emplace_back("--quiet");    // otherwise help is displayed and exited

    SECTION("--artificial-sequence-size-factor") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --artificial-sequence-size-factor 2 --dynamic-artificial-sequences");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.artificialSequenceSizeFactor() == 2);
    }

    SECTION("--dynamic-artificial-sequences") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --dynamic-artificial-sequences");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.dynamicArtificialSequences());
    }

    SECTION("--fast") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --fast");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.fast());
        REQUIRE(config.fastBatchsize() == 0);
    }

    SECTION("--fast-batchsize") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --fast --fast-batchsize 1");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.fastBatchsize() == 1);
        argvStr.pop_back(); argvStr.emplace_back("-1");  // invalid
        REQUIRE(throwsOnConstruction(argvStr));
        argvStr.pop_back(); argvStr.emplace_back("0.1");  // invalid
        REQUIRE(throwsOnConstruction(argvStr));
    }

    SECTION("--input") {
        addParameterFromString(argvStr, "--weight 8 --input");
        argvStr.emplace_back("gen1.fasta");
        REQUIRE(throwsOnConstruction(argvStr));
        argvStr.emplace_back("gen2.fasta");
        REQUIRE(!throwsOnConstruction(argvStr));
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.inputFiles() == std::vector<std::string>{"gen1.fasta","gen2.fasta"});
        REQUIRE(config.genome1() == "gen1");
        REQUIRE(config.genome2() == "gen2");
    }

    SECTION("--genome1/2") {
        addParameterFromString(argvStr, "--weight 8 --input gen1.fasta gen2.fasta gen3.fasta");
        REQUIRE(throwsOnConstruction(argvStr));
        addParameterFromString(argvStr, "--genome1 gen2");
        REQUIRE(throwsOnConstruction(argvStr));
        addParameterFromString(argvStr, "--genome2 gen3");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.genome1() == "gen2");
        REQUIRE(config.genome2() == "gen3");
    }

    SECTION("--weight-fraction") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --weight-fraction 0.4");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 20);
        REQUIRE(config.weight() == 8);
        REQUIRE(config.seedSetSize() == 1);
        addParameterFromString(argvStr, "--span 10");   // beats fraction
        config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 10);
        REQUIRE(config.weight() == 8);
        REQUIRE(config.seedSetSize() == 1);
        addParameterFromString(argvStr, "--optimal-seed");   // beats fraction and span
        config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 14);
        REQUIRE(config.weight() == 8);
        REQUIRE(config.seedSetSize() == 1);
        addParameterFromString(argvStr, "--masks 1011 1101");   // beats everything
        config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 4);
        REQUIRE(config.weight() == 3);
        REQUIRE(config.seedSetSize() == 2);
    }

    SECTION("--masks") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --masks 10101");
        REQUIRE(!throwsOnConstruction(argvStr)); // valid mask
        argvStr.emplace_back("11011");
        REQUIRE(throwsOnConstruction(argvStr)); // second mask different weight
        argvStr.pop_back(); argvStr.emplace_back("1101");
        REQUIRE(!throwsOnConstruction(argvStr)); // second mask valid
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.masks() == std::vector<std::string>{"10101","1101"});
        REQUIRE(config.span() == 5); // highest span
        REQUIRE(config.weight() == 3);
    }

    SECTION("--match-limit") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.matchLimit() == 10);    // default
        addParameterFromString(argvStr, "--match-limit 20");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.matchLimit() == 20);
        argvStr.pop_back(); argvStr.emplace_back("0");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.matchLimit() == ULLONG_MAX);  // no limit when 0
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }

    SECTION("--match-limit-discard-exceeding") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --match-limit 10");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(!config.matchLimitDiscardSeeds());    // default
        addParameterFromString(argvStr, "--match-limit-discard-exceeding");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.matchLimitDiscardSeeds());
    }

    SECTION("--weight") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 8);
        REQUIRE(config.weight() == 8);
        argvStr.pop_back(); argvStr.emplace_back("0");  // invalid
        REQUIRE(throwsOnConstruction(argvStr));
        argvStr.emplace_back("-1");                     // invalid
        REQUIRE(throwsOnConstruction(argvStr));
    }

    SECTION("--optimal-seed") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --optimal-seed --weight 33");
        REQUIRE(throwsOnConstruction(argvStr)); // no optimal for weight 33
        argvStr.pop_back(); argvStr.emplace_back("8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 14);
        REQUIRE(config.weight() == 8);
        REQUIRE(config.seedSetSize() == 1);
        addParameterFromString(argvStr, "--masks 1011 1101");   // beats optimals-seed
        config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 4);
        REQUIRE(config.weight() == 3);
        REQUIRE(config.seedSetSize() == 2);
    }

    SECTION("--output") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.output() == "seedFindingOutput.json");   // default
        std::string option{"--output "};
        std::string outTestfile{"outputTestfile.json"};
        option += outTestfile;
        auto alreadyExists = fs::exists(outTestfile);
        addParameterFromString(argvStr, option);
        config = configurationFromParameters(argvStr);
        REQUIRE(config.output() == outTestfile);
        if (!alreadyExists) {
            fs::remove(outTestfile);    // did not exist, got created, so delete again
        }
    }

    SECTION("--output-artificial-sequences") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.outputArtificialSequences() == "");   // default
        std::string option{"--output-artificial-sequences "};
        std::string outTestfile{"outputArtificialSequencesTestfile.fa"};
        option += outTestfile;
        auto alreadyExists = fs::exists(outTestfile);
        addParameterFromString(argvStr, option);
        config = configurationFromParameters(argvStr);
        REQUIRE(config.outputArtificialSequences() == outTestfile);
        if (!alreadyExists) {
            fs::remove(outTestfile);    // did not exist, got created, so delete again
        }
    }

    SECTION("--output-run-information") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.outputRunInformation() == "");   // default
        std::string option{"--output-run-information "};
        std::string outTestfile{"outputRunInformationTestfile.json"};
        option += outTestfile;
        auto alreadyExists = fs::exists(outTestfile);
        addParameterFromString(argvStr, option);
        config = configurationFromParameters(argvStr);
        REQUIRE(config.outputRunInformation() == outTestfile);
        if (!alreadyExists) {
            fs::remove(outTestfile);    // did not exist, got created, so delete again
        }
    }

    SECTION("--p") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        auto autoThreads = std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 1;
        REQUIRE(config.nThreads() == autoThreads);  // default
        addParameterFromString(argvStr, "--p 0");
        REQUIRE(throwsOnConstruction);  // explicit 0 not allowed
        argvStr.pop_back(); argvStr.emplace_back("1");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.nThreads() == 1);
    }

    SECTION("--quiet") {
        argvStr.pop_back();
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(!config.quiet());
        addParameterFromString(argvStr, "--quiet");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.quiet());
    }

    SECTION("--seed-set-size") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.seedSetSize() == 1);
        addParameterFromString(argvStr, "--seed-set-size 2");
        REQUIRE(throwsOnConstruction(argvStr)); // too many for weight == span
        argvStr.pop_back(); argvStr.pop_back();
        addParameterFromString(argvStr, "--span 10 --seed-set-size 2");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.seedSetSize() == 2);
        argvStr.pop_back(); argvStr.emplace_back("0");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }

    SECTION("--span") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 8);    // default
        REQUIRE(config.weight() == 8);
        REQUIRE(config.seedSetSize() == 1);
        addParameterFromString(argvStr, "--span 10");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 10);
        REQUIRE(config.weight() == 8);
        REQUIRE(config.seedSetSize() == 1);
        addParameterFromString(argvStr, "--optimal-seed");   // beats span
        config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 14);
        REQUIRE(config.weight() == 8);
        REQUIRE(config.seedSetSize() == 1);
        addParameterFromString(argvStr, "--masks 1011 1101");   // beats everything
        config = configurationFromParameters(argvStr);
        REQUIRE(config.span() == 4);
        REQUIRE(config.weight() == 3);
        REQUIRE(config.seedSetSize() == 2);
    }

    SECTION("--diagonal-filtering") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(!config.performDiagonalFiltering()); // default
        addParameterFromString(argvStr, "--diagonal-filtering");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.performDiagonalFiltering());
    }

    SECTION("--diagonal-threshold") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --diagonal-filtering");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.performDiagonalFiltering());
        REQUIRE(config.diagonalThreshold() == 2);    // default
        addParameterFromString(argvStr, "--diagonal-threshold 5");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.diagonalThreshold() == 5);
        REQUIRE(config.performDiagonalFiltering());
        argvStr.pop_back(); argvStr.emplace_back("0.5");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }

    SECTION("--local-search-area") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --diagonal-filtering");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.localAreaLength() == 1000);    // default
        addParameterFromString(argvStr, "--local-search-area 100");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.localAreaLength() == 100);
        argvStr.pop_back(); argvStr.emplace_back("0");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }

    SECTION("--allow-overlap") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --diagonal-filtering");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(!config.allowOverlap());    // default
        addParameterFromString(argvStr, "--allow-overlap");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.allowOverlap());
    }

    SECTION("--min-match-distance") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --diagonal-filtering");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.minMatchDistance() == 0);    // default
        addParameterFromString(argvStr, "--min-match-distance 5");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.minMatchDistance() == 5);
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr));
        argvStr.pop_back();
        addParameterFromString(argvStr, "5 --allow-overlap");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.minMatchDistance() == 0);
    }

    SECTION("--geometric-hashing") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(!config.performGeometricHashing());    // default
        addParameterFromString(argvStr, "--geometric-hashing");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.performGeometricHashing());

        REQUIRE(!config.oldCubeScore());
        addParameterFromString(argvStr, "--old-cube-score");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.oldCubeScore());
    }

    SECTION("--cube-area-cutoff") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeAreaCutoff() == 300000000);    // default
        addParameterFromString(argvStr, "--cube-area-cutoff 5");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeAreaCutoff() == 5);
        argvStr.pop_back(); argvStr.emplace_back("0");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeAreaCutoff() == 0);
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr));
    }

    SECTION("--cube-score-normalization-parameter") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreNormalizationParameter() == 300000000);    // default
        addParameterFromString(argvStr, "--cube-score-normalization-parameter 2");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreNormalizationParameter() == 2);
        argvStr.pop_back(); argvStr.emplace_back("0");
        REQUIRE(throwsOnConstruction(argvStr));
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr));
    }

    SECTION("--cube-score-parameter") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreParameter() == 500);    // default
        addParameterFromString(argvStr, "--cube-score-parameter 5");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreParameter() == 5);
        argvStr.pop_back(); argvStr.emplace_back("0");
        REQUIRE(throwsOnConstruction(argvStr));
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr));
    }

    SECTION("--cube-score-threshold") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreThreshold() == 25);    // default
        addParameterFromString(argvStr, "--cube-score-threshold 5.5");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreThreshold() == 5.5);
        argvStr.pop_back(); argvStr.emplace_back("-0.001");
        REQUIRE(throwsOnConstruction(argvStr));
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr));
    }

    SECTION("--occurrence-per-genome-max") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.occurrencePerGenomeMax() == ULLONG_MAX);    // default
        addParameterFromString(argvStr, "--occurrence-per-genome-max 10");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.occurrencePerGenomeMax() == 10);
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }

    SECTION("--occurrence-per-genome-min") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.occurrencePerGenomeMin() == 1);  // default
        addParameterFromString(argvStr, "--occurrence-per-genome-min 2");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.occurrencePerGenomeMin() == 2);
        argvStr.pop_back(); argvStr.emplace_back("0");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }

    SECTION("--tilesize") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.tileSize() == 10000);    // default
        addParameterFromString(argvStr, "--tilesize 20");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.tileSize() == 20);
        argvStr.pop_back(); argvStr.emplace_back("0");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }
}



TEST_CASE("Configuration Builder") {
    SECTION("Defaults") {
        auto configBuilder = ConfigBuilder();
        std::shared_ptr<Configuration const> config = configBuilder.makeConfig();
        REQUIRE(config->allowOverlap() == false);
        REQUIRE(config->artificialSequenceSizeFactor() == 1);
        REQUIRE(config->createAllMatches() == false);
        REQUIRE(config->cubeAreaCutoff() == 300000000);
        REQUIRE(config->cubeScoreNormalizationParameter() == 300000000);
        REQUIRE(config->cubeScoreParameter() == 500);
        REQUIRE(config->cubeScoreThreshold() == 25);
        REQUIRE(config->diagonalThreshold() == 2);
        REQUIRE(config->dynamicArtificialSequences() == false);
        REQUIRE(config->fast() == false);
        REQUIRE(config->fastBatchsize() == 0);
        REQUIRE(config->genome1() == "");
        REQUIRE(config->genome2() == "");
        REQUIRE(config->inputFiles() == std::vector<std::string>());
        REQUIRE(config->localAreaLength() == 1000);
        REQUIRE(config->masks() == std::vector<std::string>());
        REQUIRE(config->matchLimit() == 10);
        REQUIRE(config->matchLimitDiscardSeeds() == false);
        REQUIRE(config->quiet() == false);
        REQUIRE(config->nThreads() == std::thread::hardware_concurrency());
        REQUIRE(config->occurrencePerGenomeMax() == ULLONG_MAX);
        REQUIRE(config->occurrencePerGenomeMin() == 1);
        REQUIRE(config->oldCubeScore() == false);
        REQUIRE(config->optimalSeed() == false);
        REQUIRE(config->output() == "");
        REQUIRE(config->outputArtificialSequences() == "");
        REQUIRE(config->outputRunInformation() == "");
        REQUIRE(config->performDiagonalFiltering() == false);
        REQUIRE(config->performGeometricHashing() == false);
        REQUIRE(config->seedSetSize() == 1);
        REQUIRE(config->span() == 5);
        REQUIRE(config->tileSize() == 10000);
        REQUIRE(config->weight() == 5);
    }

    SECTION("Change defaults") {
        std::shared_ptr<Configuration const> config = customConfiguration(
                    AllowOverlap(true),
                    ArtificialSequenceSizeFactor(2),
                    CreateAllMatches(true),
                    CubeAreaCutoff(5),
                    CubeScoreNormalizationParameter(1),
                    CubeScoreParameter(3),
                    CubeScoreThreshold(12),
                    DiagonalThreshold(5),
                    DynamicArtificialSequences(true),
                    Fast(true),
                    FastBatchsize(32),
                    Genome1("foo"),
                    Genome2("bar"),
                    InputFiles({"file1", "file2"}),
                    LocalSearchAreaLength(50),
                    Masks({"10101", "101101"}),
                    MatchLimit(101),
                    MatchLimitDiscardSeeds(true),
                    MinMatchDistance(3),
                    NoProgressbar(true),
                    NThreads(2),
                    OccurrencePerGenomeMax(3),
                    OccurrencePerGenomeMin(2),
                    OldCubeScore(true),
                    OptimalSeed(true),
                    Output("foo.json"),
                    OutputArtificialSequences("bar.fa"),
                    OutputRunInformation("foobar.json"),
                    PerformDiagonalFiltering(true),
                    PerformGeometricHashing(true),
                    SeedSetSize(3),
                    Span(10),
                    TileSize(500),
                    Weight(9)
                    );
        REQUIRE(config->allowOverlap() == true);
        REQUIRE(config->artificialSequenceSizeFactor() == 2);
        REQUIRE(config->createAllMatches() == true);
        REQUIRE(config->cubeAreaCutoff() == 5);
        REQUIRE(config->cubeScoreNormalizationParameter() == 1);
        REQUIRE(config->cubeScoreParameter() == 3);
        REQUIRE(config->cubeScoreThreshold() == 12);
        REQUIRE(config->diagonalThreshold() == 5);
        REQUIRE(config->dynamicArtificialSequences() == true);
        REQUIRE(config->fast() == true);
        REQUIRE(config->fastBatchsize() == 32);
        REQUIRE(config->genome1() == "foo");
        REQUIRE(config->genome2() == "bar");
        REQUIRE(config->inputFiles() == std::vector<std::string>{"file1", "file2"});
        REQUIRE(config->localAreaLength() == 50);
        REQUIRE(config->masks() == std::vector<std::string>{"10101", "101101"});
        REQUIRE(config->matchLimit() == 101);
        REQUIRE(config->matchLimitDiscardSeeds() == true);
        REQUIRE(config->quiet() == true);
        REQUIRE(config->nThreads() == 2);
        REQUIRE(config->occurrencePerGenomeMax() == 3);
        REQUIRE(config->occurrencePerGenomeMin() == 2);
        REQUIRE(config->oldCubeScore() == true);
        REQUIRE(config->optimalSeed() == true);
        REQUIRE(config->output() == "foo.json");
        REQUIRE(config->outputArtificialSequences() == "bar.fa");
        REQUIRE(config->outputRunInformation() == "foobar.json");
        REQUIRE(config->performDiagonalFiltering() == true);
        REQUIRE(config->performGeometricHashing() == true);
        REQUIRE(config->seedSetSize() == 3);
        REQUIRE(config->span() == 10);
        REQUIRE(config->tileSize() == 500);
        REQUIRE(config->weight() == 9);
    }
}
