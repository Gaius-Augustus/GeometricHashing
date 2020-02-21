#include <cstdlib>
#include <limits>
#include <memory>
#include <vector>

#include "catch2/catch.hpp"
#include "../Cube.h"
#include "../Cubeset.h"
#include "../ExtractSeeds.h"
#include "../Link.h"
#include "../SeedMapContiguous.h"
#include "ConfigurationGenerator.h"
#include "GeometricHashingTestdata.h"
#include "LoadTestdata.h"

/* ~~~~~~~~~~~~~~~~~~~~
 * TESTING Tiledistance
 * ~~~~~~~~~~~~~~~~~~~~ */
TEST_CASE("Test Tiledistance") {
    // LSB needed for sign, so only allow numbers [-0x7fffffffff, 0x7fffffffff]
    auto td1 = Tiledistance(0,0,-5,true);
    auto td2 = Tiledistance(1,1,0,true);
    auto td3 = Tiledistance(2,2,5,false);
    REQUIRE(td1.distance() == -5);
    REQUIRE(td1.genome() == 0);
    REQUIRE(td1.reverse() == true);
    REQUIRE(td1.sequence() == 0);
    REQUIRE(td2.distance() == 0);
    REQUIRE(td2.genome() == 1);
    REQUIRE(td2.reverse() == true);
    REQUIRE(td2.sequence() == 1);
    REQUIRE(td3.distance() == 5);
    REQUIRE(td3.genome() == 2);
    REQUIRE(td3.reverse() == false);
    REQUIRE(td3.sequence() == 2);

    auto thrown = false;
    try {
        Tiledistance td(0,0,0x8fffffffff,false);
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        Tiledistance td(0,0,-0x8fffffffff,false);
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);

    auto td4 = Tiledistance(0,0,0x7fffffffff,false);
    auto td5 = Tiledistance(0,0,(0x7fffffffff * -1),false);
    REQUIRE(td4.distance() == 0x7fffffffff);
    REQUIRE(td5.distance() == -0x7fffffffff);

    std::set<Tiledistance> tdset;
    tdset.emplace(0,0,-5,false);
    tdset.emplace(0,0,-5,false);
    REQUIRE(tdset.size() == 1);
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~
 * TESTING Cube and Cubeset
 * ~~~~~~~~~~~~~~~~~~~~~~~~ */
TEST_CASE("Test Cube") {
    // Cube now has same dimension as Link, 0-th tiledistance is reference
    // Also, Tiledistance is only an alias for KmerOccurrence
    auto link = std::make_shared<Link>();
    link->insertOccurrence(0, 0, 10, true, "ACGT"); // ref
    link->insertOccurrence(1, 2, 30, true, "ACGT"); // same strand
    link->insertOccurrence(2, 1, 5, true, "ACGT");  // same strand, different offset
    link->insertOccurrence(3, 4, 20, false, "ACGT"); // diff. strand

    Cube cube(link, 1);
    REQUIRE(cube.dimensionality() == static_cast<uint16_t>(4));
    REQUIRE(cube.tiledistance().size() == static_cast<uint16_t>(4));

    auto& td0 = cube.tiledistance(0);
    auto& td1 = cube.tiledistance(1);
    auto& td2 = cube.tiledistance(2);
    auto& td3 = cube.tiledistance(3);

    REQUIRE(td0.distance() == static_cast<long long int>(0));
    REQUIRE(td0.genome() == static_cast<uint16_t>(0));
    REQUIRE(td0.reverse() == true);
    REQUIRE(td0.sequence() == static_cast<uint16_t>(0));

    REQUIRE(td1.distance() == static_cast<long long int>(20));
    REQUIRE(td1.genome() == static_cast<uint16_t>(1));
    REQUIRE(td1.reverse() == true);
    REQUIRE(td1.sequence() == static_cast<uint16_t>(2));

    REQUIRE(td2.distance() == static_cast<long long int>(-5));
    REQUIRE(td2.genome() == static_cast<uint16_t>(2));
    REQUIRE(td2.reverse() == true);
    REQUIRE(td2.sequence() == static_cast<uint16_t>(1));

    REQUIRE(td3.distance() == static_cast<long long int>(10));
    REQUIRE(td3.genome() == static_cast<uint16_t>(3));
    REQUIRE(td3.reverse() == false);
    REQUIRE(td3.sequence() == static_cast<uint16_t>(4));

    std::set<std::shared_ptr<Cube const>, CubePtrLess> cubeSet;
    std::unordered_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> cubeUnorderedSet;

    // test hashing / set behaviour
    auto l1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0,0,100,false,""),
                                                                 KmerOccurrence(1,1,200,false,""),
                                                                 KmerOccurrence(2,2,50,false,"")});
    auto l3 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0,0,200,false,""),
                                                                 KmerOccurrence(1,1,300,false,""),
                                                                 KmerOccurrence(2,2,150,false,"")});
    auto cube1 = std::make_shared<Cube const>(l1, 50);
    auto cube2 = std::make_shared<Cube const>(l1, 50);
    auto cube3 = std::make_shared<Cube const>(l3, 50);
    cubeSet.emplace(cube1);
    cubeSet.emplace(cube2);
    cubeSet.emplace(cube3);
    REQUIRE(cubeSet.size() == 1);
    REQUIRE(cubeSet.find(cube1) != cubeSet.end());
    REQUIRE(cubeSet.find(cube2) != cubeSet.end());
    REQUIRE(cubeSet.find(cube3) != cubeSet.end());
    cubeUnorderedSet.emplace(cube1);
    cubeUnorderedSet.emplace(cube2);
    cubeUnorderedSet.emplace(cube3);
    REQUIRE(cubeUnorderedSet.size() == 1);
    REQUIRE(cubeUnorderedSet.find(cube1) != cubeUnorderedSet.end());
    REQUIRE(cubeUnorderedSet.find(cube2) != cubeUnorderedSet.end());
    REQUIRE(cubeUnorderedSet.find(cube3) != cubeUnorderedSet.end());
}



TEST_CASE("Test Cubeset") {
    auto data = LoadTestdata(LoadTestdata::geometricHashingData);
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto fastaCollection = data.fastaCollection();

    auto p = (std::thread::hardware_concurrency() > 0)
            ? std::thread::hardware_concurrency()
            : 1;
    if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }

    //auto config = generateConfiguration("--weight 5 --cube-score-parameter 0 --cube-score-threshold 0 ");
    auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                  ArtificialSequenceSizeFactor(1),
                                                  CreateAllMatches(false),
                                                  CubeScoreMu(0),
                                                  CubeScoreThreshold(ULLONG_MAX),
                                                  DiagonalThreshold(0),
                                                  DynamicArtificialSequences(false),
                                                  Genome1(data.genome0()),
                                                  Genome2(data.genome1()),
                                                  InputFiles(data.inputFiles()),
                                                  LocalSearchAreaLength(0),
                                                  Masks(Configuration::MasksType()),
                                                  MatchLimit(ULLONG_MAX),
                                                  MatchLimitDiscardSeeds(false),
                                                  MinMatchDistance(0),
                                                  NoProgressbar(false),
                                                  NThreads(p),
                                                  OccurrencePerGenomeMax(1),
                                                  OccurrencePerGenomeMin(1),
                                                  OptimalSeed(false),
                                                  Output(""),
                                                  OutputArtificialSequences(""),
                                                  OutputRunInformation(""),
                                                  PerformDiagonalFiltering(false),
                                                  PerformGeometricHashing(true),
                                                  SeedSetSize(1),
                                                  Span(5),
                                                  TileSize(100),
                                                  Weight(5));

    auto kmerMap = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort,
                                                      TwoBitKmerDataShort>>(config,
                                                                            data.idMap(),
                                                                            std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                                                                             SpacedSeedMaskCollection::Span(config->span()),
                                                                                                                             SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize())));

    ExtractSeeds<TwoBitKmerDataShort,
                 TwoBitKmerDataShort>(fastaCollection,
                                      kmerMap,
                                      p, false);

    auto linkset = std::make_shared<Linkset>(data.idMap(), ULLONG_MAX, false, 1, 1);
    for (auto&& elem : kmerMap->seedMap()) {
        linkset->createLinks(elem.second.at(0));
    }

    auto idMap = *(linkset->idMapping()); //constIDMapRef; // true copy so that original mapping remains unchanged

    Cubeset observedCubeset(linkset, config);
    auto observedMap = observedCubeset.cubeMap();

    ManualSets manSet(idMap);
    auto& expectedCubes = manSet.expectedCubes();
    auto& expectedScores = manSet.expectedScore();    

    SECTION("Check if all/only observed are expected") {
        for (auto&& obs : observedMap) {
            REQUIRE(expectedCubes.find(obs.first) != expectedCubes.end());
        }
        for (auto&& expect : expectedCubes) {
            REQUIRE(observedMap.find(expect) != observedMap.end());
        }
    }

    SECTION("Check Cube scores") { // mu == 0, just the count of Links in the Cubes
        auto& scoreToCubeMap = observedCubeset.scoreToCube();
        for (auto&& obs : observedMap) {
            auto&& cube = obs.first;
            auto score = observedCubeset.computeCubeScore(cube);
            REQUIRE(score == expectedScores.at(cube));  // correct score
            REQUIRE(scoreToCubeMap.find(score) != scoreToCubeMap.end());    // score included in map
            REQUIRE(scoreToCubeMap.at(score).find(cube) != scoreToCubeMap.at(score).end());  // cupe associated with score
        }
    }
}

TEST_CASE("Test Cube Scoring") {
    auto p = (std::thread::hardware_concurrency() > 0)
            ? std::thread::hardware_concurrency()
            : 1;
    if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }

    auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                  ArtificialSequenceSizeFactor(1),
                                                  CreateAllMatches(false),
                                                  CubeScoreMu(3),
                                                  CubeScoreThreshold(ULLONG_MAX),
                                                  DiagonalThreshold(0),
                                                  DynamicArtificialSequences(false),
                                                  Genome1("genome0"),
                                                  Genome2("genome1"),
                                                  InputFiles(std::vector<std::string>{}),
                                                  LocalSearchAreaLength(0),
                                                  Masks(Configuration::MasksType()),
                                                  MatchLimit(ULLONG_MAX),
                                                  MatchLimitDiscardSeeds(false),
                                                  MinMatchDistance(0),
                                                  NoProgressbar(false),
                                                  NThreads(p),
                                                  OccurrencePerGenomeMax(ULLONG_MAX),
                                                  OccurrencePerGenomeMin(1),
                                                  OptimalSeed(false),
                                                  Output(""),
                                                  OutputArtificialSequences(""),
                                                  OutputRunInformation(""),
                                                  PerformDiagonalFiltering(false),
                                                  PerformGeometricHashing(true),
                                                  SeedSetSize(1),
                                                  Span(5),
                                                  TileSize(10),
                                                  Weight(5));

    auto idMap = std::make_shared<IdentifierMapping>("genome0");
    idMap->queryGenomeID("genome1");
    idMap->querySequenceID("sequence0", "genome0");
    idMap->querySequenceID("sequence1", "genome1");

    SECTION("2D Case") {
        auto linkset = std::make_shared<Linkset>(idMap,
                                                 config->matchLimit(),
                                                 config->matchLimitDiscardSeeds(),
                                                 config->occurrencePerGenomeMax(),
                                                 config->occurrencePerGenomeMin());

        auto link0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-0
                                                                        KmerOccurrence(1, 1, 1, false, "ACGT")});
        linkset->addLink(link0);
        auto cube0 = std::make_shared<Cube>(link0, config->tileSize()); // score 1, single Link

        auto link10_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-1
                                                                           KmerOccurrence(1, 1, 11, false, "ACGT")});
        auto link10_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 10, false, "ACGT"),
                                                                           KmerOccurrence(1, 1, 21, false, "ACGT")});
        linkset->addLink(link10_0);
        linkset->addLink(link10_1);
        auto cube1 = std::make_shared<Cube>(link10_0, config->tileSize()); // score 8, two Links on same diagonal

        auto link20_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-2
                                                                           KmerOccurrence(1, 1, 21, false, "ACGT")});
        auto link20_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),
                                                                           KmerOccurrence(1, 1, 21, false, "ACGT")});
        linkset->addLink(link20_0);
        linkset->addLink(link20_1);
        auto cube2 = std::make_shared<Cube>(link20_0, config->tileSize()); // score 8, single Link with count 2

        auto link30_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-3
                                                                           KmerOccurrence(1, 1, 31, false, "ACGT")});
        auto link30_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 1, false, "ACGT"),        // diag +1
                                                                           KmerOccurrence(1, 1, 33, false, "ACGT")});
        auto link30_2 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 2, false, "ACGT"),
                                                                           KmerOccurrence(1, 1, 34, false, "ACGT")});
        auto link30_3 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 3, false, "ACGT"),        // diag +3
                                                                           KmerOccurrence(1, 1, 37, false, "ACGT")});
        linkset->addLink(link30_0);
        linkset->addLink(link30_1);
        linkset->addLink(link30_2);
        linkset->addLink(link30_3);
        auto cube3 = std::make_shared<Cube>(link30_0, config->tileSize()); // score 12.5, one delta 1, two Links on same diagonal, one delta 2

        auto link40_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-4
                                                                           KmerOccurrence(1, 1, 41, false, "ACGT")});
        auto link40_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 1, false, "ACGT"),        // diag +1
                                                                           KmerOccurrence(1, 1, 43, false, "ACGT")});
        auto link40_2 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 2, false, "ACGT"),
                                                                           KmerOccurrence(1, 1, 44, false, "ACGT")});
        auto link40_3 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 2, false, "ACGT"),
                                                                           KmerOccurrence(1, 1, 44, false, "ACGT")});
        auto link40_4 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 3, false, "ACGT"),        // diag +3
                                                                           KmerOccurrence(1, 1, 47, false, "ACGT")});
        linkset->addLink(link40_0);
        linkset->addLink(link40_1);
        linkset->addLink(link40_2);
        linkset->addLink(link40_3);
        linkset->addLink(link40_4);
        auto cube4 = std::make_shared<Cube>(link40_0, config->tileSize()); // score 16.5, one delta 1, two Links on same diagonal, one of them with count 2, one delta 2

        auto link50_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-5, diag -1
                                                                           KmerOccurrence(1, 1, 51, false, "ACGT")});
        auto link50_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 3, false, "ACGT"),        // diag +1 --> now goes into cube 4 as diag +8
                                                                           KmerOccurrence(1, 1, 52, false, "ACGT")});
        linkset->addLink(link50_0);
        linkset->addLink(link50_1);
        auto cube5 = std::make_shared<Cube>(link50_0, config->tileSize()); // score 4, diag difference is 2 for both Links --> should be 1

        auto cubeset = Cubeset(linkset, config);
        REQUIRE(cubeset.computeCubeScore(cube0) == 1.);
        REQUIRE(cubeset.computeCubeScore(cube1) == 8.);
        REQUIRE(cubeset.computeCubeScore(cube2) == 8.);
        REQUIRE(cubeset.computeCubeScore(cube3) == 12.5);
        REQUIRE(cubeset.computeCubeScore(cube4) == 18.);//16.5);
        REQUIRE(cubeset.computeCubeScore(cube5) == 1.);//4);
    }

    SECTION("Higher Dimensional Case") {
        idMap->queryGenomeID("genome2");
        idMap->queryGenomeID("genome3");
        idMap->queryGenomeID("genome4");
        idMap->querySequenceID("sequence2", "genome2");
        idMap->querySequenceID("sequence3", "genome3");
        idMap->querySequenceID("sequence4", "genome4");
        auto linkset = std::make_shared<Linkset>(idMap,
                                                 config->matchLimit(),
                                                 config->matchLimitDiscardSeeds(),
                                                 config->occurrencePerGenomeMax(),
                                                 config->occurrencePerGenomeMin());

        auto link0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-0-0-0-0
                                                                        KmerOccurrence(1, 1, 0, false, "ACGT"),
                                                                        KmerOccurrence(2, 2, 0, false, "ACGT"),
                                                                        KmerOccurrence(3, 3, 0, false, "ACGT"),
                                                                        KmerOccurrence(4, 4, 0, false, "ACGT")});
        linkset->addLink(link0);
        auto cube0 = std::make_shared<Cube>(link0, config->tileSize()); // score 1, single Link

        auto link1_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-1-0-0-0
                                                                          KmerOccurrence(1, 1, 10, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 0, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 0, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 0, false, "ACGT")});
        auto link1_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 10, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 20, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 10, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 10, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 10, false, "ACGT")});
        linkset->addLink(link1_0);
        linkset->addLink(link1_1);
        auto cube1 = std::make_shared<Cube>(link1_0, config->tileSize()); // score 8, two Links on same diagonal

        auto link2_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 10, false, "ACGT"),    // cube 0-1-1-0-0
                                                                          KmerOccurrence(1, 1, 20, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 20, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 10, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 10, false, "ACGT")});
        auto link2_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 10, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 20, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 20, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 10, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 10, false, "ACGT")});
        linkset->addLink(link2_0);
        linkset->addLink(link2_1);
        auto cube2 = std::make_shared<Cube>(link2_0, config->tileSize()); // score 8, one Link with count 2

        auto link3_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-1-1-1-0
                                                                          KmerOccurrence(1, 1, 10, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 10, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 10, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 0, false, "ACGT")});
        auto link3_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 10, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 21, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 21, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 21, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 10, false, "ACGT")});
        auto link3_2 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 11, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 22, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 22, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 22, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 11, false, "ACGT")});
        auto link3_3 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 20, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 33, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 33, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 33, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 20, false, "ACGT")});
        linkset->addLink(link3_0);
        linkset->addLink(link3_1);
        linkset->addLink(link3_2);
        linkset->addLink(link3_3);
        auto mu = config->cubeScoreMu();
        auto score3 = 1. + (mu / (1. + std::sqrt(3)))
                      + 1. + (mu / (1. + std::sqrt(0)))
                      + 1. + (mu / (1. + std::sqrt(0)))
                      + 1. + (mu / (1. + std::sqrt(12)));
        auto cube3 = std::make_shared<Cube>(link3_0, config->tileSize()); // score ?, one Link with delta ?, two Links on same diagonal, one Link with delta ?

        auto link4_0 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-1-1-1-1
                                                                          KmerOccurrence(1, 1, 13, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 13, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 13, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 13, false, "ACGT")});
        auto link4_1 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 10, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 22, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 22, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 22, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 22, false, "ACGT")});
        auto link4_2 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 11, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 23, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 23, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 23, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 23, false, "ACGT")});
        auto link4_3 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 11, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 23, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 23, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 23, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 23, false, "ACGT")});
        auto link4_4 = std::make_shared<Link>(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 20, false, "ACGT"),
                                                                          KmerOccurrence(1, 1, 30, false, "ACGT"),
                                                                          KmerOccurrence(2, 2, 30, false, "ACGT"),
                                                                          KmerOccurrence(3, 3, 30, false, "ACGT"),
                                                                          KmerOccurrence(4, 4, 30, false, "ACGT")});
        linkset->addLink(link4_0);
        linkset->addLink(link4_1);
        linkset->addLink(link4_2);
        linkset->addLink(link4_3);
        linkset->addLink(link4_4);
        auto score4 = 1. + (mu / (1. + std::sqrt(4)))
                      + 1. + (mu / (1. + std::sqrt(0)))
                      + (2. * (mu + 1.))
                      + 1. + (mu / (1. + std::sqrt(16)));
        auto cube4 = std::make_shared<Cube>(link4_0, config->tileSize()); // score ?, one Link with delta ?, two Links on same diagonal, one Link with delta ?

        auto cubeset = Cubeset(linkset, config);
        REQUIRE(cubeset.computeCubeScore(cube0) == 1.);
        REQUIRE(cubeset.computeCubeScore(cube1) == 8.);
        REQUIRE(cubeset.computeCubeScore(cube2) == 8.);
        REQUIRE(cubeset.computeCubeScore(cube3) == score3);
        REQUIRE(cubeset.computeCubeScore(cube4) == score4);
    }
}

TEST_CASE("Test Link Grouping") {
    auto link1 = std::make_shared<Link>();
    link1->insertOccurrence(0, 0, 0, false, "ACGT");
    link1->insertOccurrence(1, 2, 0, false, "ACGT");
    link1->insertOccurrence(2, 4, 10, false, "ACGT");

    auto link2 = std::make_shared<Link>();  // same cube
    link2->insertOccurrence(0, 0, 10, false, "ACGT");
    link2->insertOccurrence(1, 2, 10, false, "ACGT");
    link2->insertOccurrence(2, 4, 20, false, "ACGT");

    auto link3 = std::make_shared<Link>();  // different cube
    link3->insertOccurrence(0, 0, 0, false, "ACGT");
    link3->insertOccurrence(1, 1, 0, false, "ACGT");
    link3->insertOccurrence(2, 4, 10, false, "ACGT");

    auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                  ArtificialSequenceSizeFactor(1),
                                                  CreateAllMatches(false),
                                                  CubeScoreMu(3),
                                                  CubeScoreThreshold(1),
                                                  DiagonalThreshold(0),
                                                  DynamicArtificialSequences(false),
                                                  Genome1(""),
                                                  Genome2(""),
                                                  InputFiles({""}),
                                                  LocalSearchAreaLength(0),
                                                  Masks(Configuration::MasksType()),
                                                  MatchLimit(ULLONG_MAX),
                                                  MatchLimitDiscardSeeds(false),
                                                  MinMatchDistance(0),
                                                  NoProgressbar(false),
                                                  NThreads(1),
                                                  OccurrencePerGenomeMax(ULLONG_MAX),
                                                  OccurrencePerGenomeMin(0),
                                                  OptimalSeed(false),
                                                  Output(""),
                                                  OutputArtificialSequences(""),
                                                  OutputRunInformation(""),
                                                  PerformDiagonalFiltering(false),
                                                  PerformGeometricHashing(true),
                                                  SeedSetSize(1),
                                                  Span(5),
                                                  TileSize(1),
                                                  Weight(5));

    auto idMapping = std::make_shared<IdentifierMapping>("ref");
    auto id = idMapping->queryGenomeID("g1");
    REQUIRE(id == 1);
    id = idMapping->queryGenomeID("g2");
    REQUIRE(id == 2);
    id = idMapping->querySequenceID("s0", "ref");
    REQUIRE(id == 0);
    id = idMapping->querySequenceID("s1", "g1");
    REQUIRE(id == 1);
    id = idMapping->querySequenceID("s2", "g2");
    REQUIRE(id == 2);
    id = idMapping->querySequenceID("s3", "ref");
    REQUIRE(id == 3);
    id = idMapping->querySequenceID("s4", "ref");
    REQUIRE(id == 4);

    auto linkset = std::make_shared<Linkset>(idMapping,
                                             config->matchLimit(),
                                             config->matchLimitDiscardSeeds(),
                                             config->occurrencePerGenomeMax(),
                                             config->occurrencePerGenomeMin());
    linkset->addLink(link1);
    linkset->addLink(link1);
    linkset->addLink(link2);
    linkset->addLink(link3);
    REQUIRE(linkset->linkset().find(link1) != linkset->linkset().end());
    REQUIRE(linkset->linkset().find(link2) != linkset->linkset().end());
    REQUIRE(linkset->linkset().find(link3) != linkset->linkset().end());
    REQUIRE(linkset->linkset().at(link1) == 2); // Link count
    REQUIRE(linkset->linkset().at(link2) == 1);
    REQUIRE(linkset->linkset().at(link3) == 1);



    Cubeset cubeset(linkset, config);
    auto cube1 = std::make_shared<Cube>(link1, config->tileSize());
    auto cube2 = std::make_shared<Cube>(link2, config->tileSize());
    auto cube3 = std::make_shared<Cube>(link3, config->tileSize());

    CubePtrEqual ceq;
    REQUIRE(ceq(cube1, cube1));
    REQUIRE(ceq(cube1, cube2));
    REQUIRE(!ceq(cube1, cube3));
    REQUIRE(ceq(cube2, cube2));
    REQUIRE(!ceq(cube2, cube3));
    REQUIRE(ceq(cube3, cube3));

    REQUIRE(cubeset.cubeMap().size() == 2);
    REQUIRE(cubeset.links(cube1).size() == 2);
    REQUIRE(cubeset.links(cube2).size() == 2);
    REQUIRE(cubeset.links(cube3).size() == 1);

    REQUIRE(cubeset.links(cube1) == cubeset.links(cube2));
    REQUIRE(cubeset.links(cube1) != cubeset.links(cube3));
    REQUIRE(cubeset.links(cube1).find(link1) != cubeset.links(cube1).end());
    REQUIRE(cubeset.links(cube1).find(link2) != cubeset.links(cube1).end());
    REQUIRE(cubeset.links(cube1).find(link3) == cubeset.links(cube1).end());
    REQUIRE(cubeset.links(cube2).find(link1) != cubeset.links(cube2).end());
    REQUIRE(cubeset.links(cube2).find(link2) != cubeset.links(cube2).end());
    REQUIRE(cubeset.links(cube2).find(link3) == cubeset.links(cube2).end());
    REQUIRE(cubeset.links(cube3).find(link1) == cubeset.links(cube3).end());
    REQUIRE(cubeset.links(cube3).find(link2) == cubeset.links(cube3).end());
    REQUIRE(cubeset.links(cube3).find(link3) != cubeset.links(cube3).end());
}



TEST_CASE("Test Link Diagonal Sorting") {
    auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                  ArtificialSequenceSizeFactor(1),
                                                  CreateAllMatches(false),
                                                  CubeScoreMu(3),
                                                  CubeScoreThreshold(1),
                                                  DiagonalThreshold(0),
                                                  DynamicArtificialSequences(false),
                                                  Genome1(""),
                                                  Genome2(""),
                                                  InputFiles({""}),
                                                  LocalSearchAreaLength(0),
                                                  Masks(Configuration::MasksType()),
                                                  MatchLimit(ULLONG_MAX),
                                                  MatchLimitDiscardSeeds(false),
                                                  MinMatchDistance(0),
                                                  NoProgressbar(false),
                                                  NThreads(1),
                                                  OccurrencePerGenomeMax(ULLONG_MAX),
                                                  OccurrencePerGenomeMin(0),
                                                  OptimalSeed(false),
                                                  Output(""),
                                                  OutputArtificialSequences(""),
                                                  OutputRunInformation(""),
                                                  PerformDiagonalFiltering(false),
                                                  PerformGeometricHashing(true),
                                                  SeedSetSize(1),
                                                  Span(5),
                                                  TileSize(20),
                                                  Weight(5));

    auto link1 = std::make_shared<Link>();
    link1->insertOccurrence(0, 0, 0, false, "ACGT");
    link1->insertOccurrence(1, 2, 0, false, "ACGT");
    link1->insertOccurrence(2, 4, 10, false, "ACGT");

    auto link2 = std::make_shared<Link>();
    link2->insertOccurrence(0, 0, 1, false, "ACGT");
    link2->insertOccurrence(1, 2, 1, false, "ACGT");
    link2->insertOccurrence(2, 4, 11, false, "ACGT");

    auto link3 = std::make_shared<Link>();
    link3->insertOccurrence(0, 0, 2, false, "ACGT");
    link3->insertOccurrence(1, 2, 3, false, "ACGT");    // more negative distance
    link3->insertOccurrence(2, 4, 12, false, "ACGT");

    auto link4 = std::make_shared<Link>();
    link4->insertOccurrence(0, 0, 3, false, "ACGT");
    link4->insertOccurrence(1, 2, 3, false, "ACGT");
    link4->insertOccurrence(2, 4, 13, false, "ACGT");

    auto link5 = std::make_shared<Link>();
    link5->insertOccurrence(0, 0, 4, false, "ACGT");
    link5->insertOccurrence(1, 2, 4, false, "ACGT");
    link5->insertOccurrence(2, 4, 12, false, "ACGT");   // more positive distance

    auto link6 = std::make_shared<Link>();
    link6->insertOccurrence(0, 0, 5, false, "ACGT");
    link6->insertOccurrence(1, 2, 6, false, "ACGT");    // more negative
    link6->insertOccurrence(2, 4, 16, false, "ACGT");   // also more negative

    auto expectedOrder = std::vector<std::shared_ptr<Link const>>{link6, link3, link1, link2, link4, link5};

    // create cubeset
    auto idMapping = std::make_shared<IdentifierMapping>("ref");
    //auto id = idMapping->queryGenomeID("g1");
    idMapping->queryGenomeID("g2");
    idMapping->querySequenceID("s0", "ref");
    idMapping->querySequenceID("s1", "g1");
    idMapping->querySequenceID("s2", "g2");
    idMapping->querySequenceID("s3", "ref");
    idMapping->querySequenceID("s4", "ref");

    auto linkset = std::make_shared<Linkset>(idMapping,
                                             config->matchLimit(),
                                             config->matchLimitDiscardSeeds(),
                                             config->occurrencePerGenomeMax(),
                                             config->occurrencePerGenomeMin());
    linkset->addLink(link1);
    linkset->addLink(link1);
    linkset->addLink(link2);
    linkset->addLink(link3);
    linkset->addLink(link4);
    linkset->addLink(link5);
    linkset->addLink(link6);

    Cubeset cubeset(linkset, config);
    auto cube = std::make_shared<Cube>(link1, config->tileSize());

    REQUIRE(cubeset.cubeMap().size() == 1);
    REQUIRE(cubeset.cubeMap().find(cube) != cubeset.cubeMap().end());
    auto& setOfLinks = cubeset.cubeMap().at(cube);
    auto observedOrder = std::vector<std::shared_ptr<Link const>>{setOfLinks.begin(), setOfLinks.end()};
    REQUIRE(observedOrder == expectedOrder);
}
