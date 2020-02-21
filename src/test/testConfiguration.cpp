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
        REQUIRE(config.matchLimit() == ULLONG_MAX);    // default
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

    SECTION("--diagonal-threshold") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.diagonalThreshold() == 0);    // default
        REQUIRE(!config.performDiagonalFiltering());
        addParameterFromString(argvStr, "--diagonal-threshold 2");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.diagonalThreshold() == 2);
        REQUIRE(config.performDiagonalFiltering());
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }

    SECTION("--local-search-area") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --diagonal-threshold 2");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.localAreaLength() == 100);    // default
        addParameterFromString(argvStr, "--local-search-area 1000");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.localAreaLength() == 1000);
        argvStr.pop_back(); argvStr.emplace_back("0");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid for user
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }

    SECTION("--allow-overlap") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --diagonal-threshold 2");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(!config.allowOverlap());    // default
        addParameterFromString(argvStr, "--allow-overlap");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.allowOverlap());
    }

    SECTION("--min-match-distance") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --diagonal-threshold 2");
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
    }

    SECTION("--cube-score-parameter") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreMu() == 3);    // default
        addParameterFromString(argvStr, "--cube-score-parameter 5.5");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreMu() == 5.5);
        argvStr.pop_back(); argvStr.emplace_back("-0.001");
        REQUIRE(throwsOnConstruction(argvStr));
    }

    SECTION("--cube-score-threshold") {
        addParameterFromString(argvStr, "--input gen1.fasta gen2.fasta --weight 8 --geometric-hashing");
        auto config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreThreshold() == 10);    // default
        addParameterFromString(argvStr, "--cube-score-threshold 5");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.cubeScoreThreshold() == 5);
        argvStr.pop_back(); argvStr.emplace_back("-0.001");
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
        REQUIRE(config.tileSize() == 100);    // default
        addParameterFromString(argvStr, "--tilesize 20");
        config = configurationFromParameters(argvStr);
        REQUIRE(config.tileSize() == 20);
        argvStr.pop_back(); argvStr.emplace_back("0");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
        argvStr.pop_back(); argvStr.emplace_back("-1");
        REQUIRE(throwsOnConstruction(argvStr)); // invalid
    }
}
