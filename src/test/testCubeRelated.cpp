#include <cstdlib>
#include <limits>
#include <memory>
#include <vector>

#include "catch2/catch.hpp"
#include "../Cube.h"
#include "../Cubeset.h"
#include "../Linkset.h"
#include "ConfigurationGenerator.h"
#include "GeometricHashingTestdata.h"
#include "LoadTestdata.h"

/* ~~~~~~~~~~~~~~~~~~~~
 * TESTING Tiledistance
 * ~~~~~~~~~~~~~~~~~~~~ */

TEST_CASE("Test Tiledistance") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Tiledistance" << std::endl;
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



/* ~~~~~~~~~~~~
 * TESTING Cube
 * ~~~~~~~~~~~~ */

TEST_CASE("Test Cube") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Cube" << std::endl;
    // Cube now has same dimension as Link, 0-th tiledistance is reference
    // Also, Tiledistance is only an alias for KmerOccurrence
    auto link = Link();
    link.insertOccurrence(0, 0, 10, true, "ACGT"); // ref
    link.insertOccurrence(1, 2, 30, true, "ACGT"); // same strand
    link.insertOccurrence(2, 1, 5, true, "ACGT");  // same strand, different offset
    link.insertOccurrence(3, 4, 20, false, "ACGT"); // diff. strand

    Cube cube(link, 1);
    REQUIRE(cube.dimensionality() == static_cast<uint16_t>(4));
    REQUIRE(cube.tiledistance().size() == static_cast<uint16_t>(4));

    auto& td0 = cube.tiledistance(0);
    auto& td1 = cube.tiledistance(1);
    auto& td2 = cube.tiledistance(2);
    auto& td3 = cube.tiledistance(3);
    REQUIRE(cube.tileindex(0) == td0.distance());
    REQUIRE(cube.tileindex(1) == td1.distance());
    REQUIRE(cube.tileindex(2) == td2.distance());
    REQUIRE(cube.tileindex(3) == td3.distance());

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

    cube = Cube(link, 5);
    REQUIRE(cube.tileindex(0) == 0);
    REQUIRE(cube.tileindex(1) == 4);
    REQUIRE(cube.tileindex(2) == -1);
    REQUIRE(cube.tileindex(3) == 2);



    std::set<std::shared_ptr<Cube const>, CubePtrLess> cubeSet;
    std::unordered_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> cubeUnorderedSet;

    // test hashing / set behaviour
    auto l1 = Link(std::vector<KmerOccurrence>{KmerOccurrence(0,0,100,false,""),
                                               KmerOccurrence(1,1,200,false,""),
                                               KmerOccurrence(2,2,50,false,"")}, 4);
    auto l3 = Link(std::vector<KmerOccurrence>{KmerOccurrence(0,0,200,false,""),
                                               KmerOccurrence(1,1,300,false,""),
                                               KmerOccurrence(2,2,150,false,"")}, 4);
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

    REQUIRE(cube.positionsToTile(0,10,10) == -1);
    REQUIRE(cube.positionsToTile(0,9,10) == -1);
    REQUIRE(cube.positionsToTile(0,11,10) == -2);
    REQUIRE(cube.positionsToTile(0,0,10) == 0);
    REQUIRE(cube.positionsToTile(9,0,10) == 0);
    REQUIRE(cube.positionsToTile(10,0,10) == 1);
    REQUIRE(cube.positionsToTile(10,1,10) == 0);
    REQUIRE(cube.positionsToTile(10,9,10) == 0);
    REQUIRE(cube.positionsToTile(11,10,10) == 0);
    REQUIRE(cube.positionsToTile(11,0,10) == 1);
    REQUIRE(cube.positionsToTile(19,0,10) == 1);
    REQUIRE(cube.positionsToTile(20,0,10) == 2);
}

TEST_CASE("Test Cube Geo and Chunks") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Cube Geo and Chunks" << std::endl;
    size_t tilesize = 12;
    size_t chunksize = 7;

    tsl::hopscotch_map<size_t, size_t> seqlens;
    seqlens[0] = 40;
    seqlens[1] = 47;
    seqlens[2] = 35;
    seqlens[3] = 57;

    size_t total = 1; for (auto&& elem : seqlens) { total *= elem.second; }
    tsl::hopscotch_map<Cube, std::array<size_t, 4>, CubeHash> cubeToGeo; // len, minchunk, maxchunk, linkcount = volume
    mabl3::ProgressBar pb{total};
    for (size_t i = 0; i < seqlens[0]; ++i) {
        for (size_t j = 0; j < seqlens[1]; ++j) {
            for (size_t k = 0; k < seqlens[2]; ++k) {
                for (size_t l = 0; l < seqlens[3]; ++l) {
                    Link link(std::vector<KmerOccurrence>{
                                  KmerOccurrence{0,0,i,false,""},
                                  KmerOccurrence{1,1,j,false,""},
                                  KmerOccurrence{2,2,k,false,""},
                                  KmerOccurrence{3,3,l,false,""}},1);
                    Cube cube(link, tilesize);
                    REQUIRE(static_cast<size_t>(std::floor(static_cast<double>(i+j+k+l) / static_cast<double>(chunksize))) == link.chunkID(chunksize));
                    if (cubeToGeo.find(cube) == cubeToGeo.end()) { cubeToGeo.insert({cube, std::array<size_t,4>{0,SIZE_MAX,0,0}}); }
                    if (((long long)tilesize * cube.tiledistance(1).distance() == (long long)j - (long long)i)
                            && ((long long)tilesize * cube.tiledistance(2).distance() == (long long)k - (long long)i)
                            && ((long long)tilesize * cube.tiledistance(3).distance() == (long long)l - (long long)i)) {
                        cubeToGeo[cube].at(0) += 1; // count length on cube root
                    }
                    // remember smallest/biggest chunks
                    auto chunk = link.chunkID(chunksize);
                    cubeToGeo[cube].at(1) = std::min(cubeToGeo[cube].at(1), chunk);
                    cubeToGeo[cube].at(2) = std::max(cubeToGeo[cube].at(2), chunk);
                    cubeToGeo[cube].at(3) += 1;
                    ++pb;
                }
            }
        }
    }
    for (auto&& elem : cubeToGeo) {
        auto geo = elem.first.geoProperties(tilesize, chunksize, seqlens);
        //std::cout << "[DEBUG] -- cube " << elem.first << ", len " << geo.length << ", min " << geo.minChunk << ", max " << geo.maxChunk << ", vol " << elem.second.at(3) << std::endl;
        REQUIRE(elem.second.at(0) == geo.length);
        if (geo.length > 0) {
            REQUIRE(elem.second.at(1) == geo.minChunk);
            REQUIRE(elem.second.at(2) == geo.maxChunk);
        } else {
            REQUIRE(elem.second.at(1) >= geo.minChunk);
            REQUIRE(elem.second.at(2) <= geo.maxChunk);
        }
    }
}



/* ~~~~~~~~~~~~~~~
 * TESTING Cubeset
 * ~~~~~~~~~~~~~~~ */

TEST_CASE("Test Cubeset") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Cubeset" << std::endl;
    auto data = LoadTestdata(LoadTestdata::geometricHashingData);
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto fastaCollection = data.fastaCollection();
    auto sequenceLengths = data.sequenceLengths();
    auto config = customConfiguration(CubeScoreParameter(3),
                                      CubeScoreThreshold(ULLONG_MAX),
                                      Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      Hasse(true),
                                      InputFiles(data.inputFiles()),
                                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                         SpacedSeedMaskCollection::Span{5},
                                                                                                         SpacedSeedMaskCollection::SeedSetSize{1})},
                                      MatchLimit{SIZE_MAX},
                                      MatchLimitDiscardSeeds{false},
                                      OccurrencePerGenomeMax{1},
                                      OccurrencePerGenomeMin{1},
                                      OldCubeScore(false),
                                      PerformGeometricHashing(true),
                                      Thinning(1),
                                      TileSize(100));
    auto kmerMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, data.idMap());
    kmerMap->extractSeeds(fastaCollection);

    auto linkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, kmerMap->idMap());
    linkset->createLinks(*kmerMap);

    auto idMap = *(linkset->idMapping()); //constIDMapRef; // true copy so that original mapping remains unchanged

    Cubeset observedCubeset(linkset, sequenceLengths);
    REQUIRE(observedCubeset.nsubtiles() == 3);
    REQUIRE(observedCubeset.numGenomes() == 5);
    auto observedMap = observedCubeset.cubeMap();

    ManualSets manSet(idMap);
    auto& expectedCubes = manSet.expectedCubes();
    // auto& expectedScores = manSet.expectedScore();

    SECTION("Check if all/only observed are expected") {
        for (auto&& obs : observedMap) {
            REQUIRE(expectedCubes.find(obs.first) != expectedCubes.end());
        }
        for (auto&& expect : expectedCubes) {
            REQUIRE(observedMap.find(expect) != observedMap.end());
        }
    }
/*
    SECTION("Area Estimation") {
        // a little bootstrapping, mainly checking if poperties are queried correctly
        auto tileSize = static_cast<double>(config->tileSize());
        SECTION("S1 >> S2") {
            //auto g1name = "species1";
            auto s1name = "species1_sequence1";
            double s1len = 111;
            //auto g2name = "species3";
            auto s2name = "species3_sequence1";
            double s2len = 6;
            for (double i = 0; i < s1len; ++i ) {
                for (double j = 0; j < s2len; ++j) {
                    auto link = Link();
                    link.insertOccurrence(idMap.queryGenomeIDConst("species1"),
                                          idMap.querySequenceID_LEGACY(s1name),
                                          i, false, "ACGT");
                    link.insertOccurrence(idMap.queryGenomeIDConst("species3"),
                                          idMap.querySequenceID_LEGACY(s2name),
                                          j, false, "ACGT");
                    auto cube = Cube(link, static_cast<size_t>(tileSize));
                    auto dist = std::floor((j-i)/tileSize);
                    auto v = std::min(s2len - (tileSize*dist), s1len);
                    auto u = std::max(-1*tileSize*dist, 0.);
                    auto A = (v - u) * tileSize;
                    REQUIRE(observedCubeset.area(cube) == static_cast<size_t>(std::max(A,0.)));
                }
            }
        }
        SECTION("S1 ~ S2") {
            auto s1name = "species1_sequence1";
            double s1len = 111;
            auto s2name = "species2_sequence1";
            double s2len = 111;
            for (double i = 0; i < s1len; ++i ) {
                for (double j = 0; j < s2len; ++j) {
                    auto link = Link();
                    link.insertOccurrence(idMap.queryGenomeIDConst("species1"),
                                          idMap.querySequenceID_LEGACY(s1name),
                                          i, false, "ACGT");
                    link.insertOccurrence(idMap.queryGenomeIDConst("species3"),
                                          idMap.querySequenceID_LEGACY(s2name),
                                          j, false, "ACGT");
                    auto cube = Cube(link, static_cast<size_t>(tileSize));
                    auto dist = std::floor((j-i)/tileSize);
                    auto v = std::min(s2len - (tileSize*dist), s1len);
                    auto u = std::max(-1*tileSize*dist, 0.);
                    auto A = (v - u) * tileSize;
                    REQUIRE(observedCubeset.area(cube) == static_cast<size_t>(std::max(A,0.)));
                }
            }
        }
        SECTION("S1 << S2") {
            auto s1name = "species4_sequence1";
            double s1len = 12;
            auto s2name = "species2_sequence1";
            double s2len = 111;
            for (double i = 0; i < s1len; ++i ) {
                for (double j = 0; j < s2len; ++j) {
                    auto link = Link();
                    link.insertOccurrence(idMap.queryGenomeIDConst("species1"),
                                          idMap.querySequenceID_LEGACY(s1name),
                                          i, false, "ACGT");
                    link.insertOccurrence(idMap.queryGenomeIDConst("species3"),
                                          idMap.querySequenceID_LEGACY(s2name),
                                          j, false, "ACGT");
                    auto cube = Cube(link, static_cast<size_t>(tileSize));
                    auto dist = std::floor((j-i)/tileSize);
                    auto v = std::min(s2len - (tileSize*dist), s1len);
                    auto u = std::max(-1*tileSize*dist, 0.);
                    auto A = (v - u) * tileSize;
                    REQUIRE(observedCubeset.area(cube) == static_cast<size_t>(std::max(A,0.)));
                }
            }
        }
    }*/
}

TEST_CASE("Test Subcube Finding") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Cube Scoring" << std::endl;
    auto config = customConfiguration(Genome1("genome0"),
                                      Genome2("genome1"),
                                      Hasse(true),
                                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                         SpacedSeedMaskCollection::Span{5},
                                                                                                         SpacedSeedMaskCollection::SeedSetSize{1})},
                                      NThreads(1));
    auto idMap = std::make_shared<IdentifierMapping>("genome0");
    idMap->queryGenomeID("genome1");
    idMap->queryGenomeID("genome2");
    idMap->queryGenomeID("genome3");
    idMap->queryGenomeID("genome4");
    idMap->querySequenceID("sequence0", "genome0");
    idMap->querySequenceID("sequence1", "genome1");
    idMap->querySequenceID("sequence2", "genome2");
    idMap->querySequenceID("sequence3", "genome3");
    idMap->querySequenceID("sequence4", "genome4");
    idMap->querySequenceID("sequence5", "genome2");
    auto seqLens = std::make_shared<tsl::hopscotch_map<size_t, size_t>>();
    (*seqLens)[0] = 1000;
    (*seqLens)[1] = 1500;
    (*seqLens)[2] = 700;
    (*seqLens)[3] = 1700;
    (*seqLens)[4] = 2000;
    (*seqLens)[5] = 1100;

    auto linkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMap);
    Cubeset::SubcubeMapType expectedSubcubes{};

    auto link = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),
                                                    KmerOccurrence(1, 1, 0, false, "ACGT"),
                                                    KmerOccurrence(2, 2, 0, false, "ACGT"),
                                                    KmerOccurrence(3, 3, 0, false, "ACGT"),
                                                    KmerOccurrence(4, 4, 0, false, "ACGT")}, 5);
    linkset->addLink(link);
    auto cube0 = std::make_shared<Cube>(*link, config->tileSize());

    link = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),
                                               KmerOccurrence(2, 2, 0, false, "ACGT"),
                                               KmerOccurrence(3, 3, 0, false, "ACGT"),
                                               KmerOccurrence(4, 4, 0, false, "ACGT")}, 5); // subcube of cube0
    linkset->addLink(link);
    auto cube1 = std::make_shared<Cube>(*link, config->tileSize());

    link = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),
                                               KmerOccurrence(1, 1, 0, false, "ACGT"),
                                               KmerOccurrence(2, 2, 0, false, "ACGT"),
                                               KmerOccurrence(3, 3, 0, false, "ACGT")}, 5); // subcube of cube0
    linkset->addLink(link);
    auto cube2 = std::make_shared<Cube>(*link, config->tileSize());

    link = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),
                                               KmerOccurrence(2, 5, 0, false, "ACGT"),
                                               KmerOccurrence(3, 3, 0, false, "ACGT"),
                                               KmerOccurrence(4, 4, 0, false, "ACGT")}, 5); // not a subcube
    linkset->addLink(link);
    auto cube3 = std::make_shared<Cube>(*link, config->tileSize());

    link = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),
                                               KmerOccurrence(1, 1, 0, false, "ACGT"),
                                               KmerOccurrence(3, 3, 0, false, "ACGT")}, 5); // subcube of cube0, cube2
    linkset->addLink(link);
    auto cube4 = std::make_shared<Cube>(*link, config->tileSize());

    link = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),
                                               KmerOccurrence(4, 4, 0, false, "ACGT")}, 5); // subcube of c0, c1, c3
    linkset->addLink(link);
    auto cube5 = std::make_shared<Cube>(*link, config->tileSize());





    expectedSubcubes[cube0].insert({cube1, cube2, cube4, cube5});
    expectedSubcubes[cube1].insert({cube5});
    expectedSubcubes[cube2].insert({cube4});
    expectedSubcubes[cube3].insert({cube5});

    REQUIRE(!(cube0->hasSubcube(*cube0)));
    REQUIRE(cube0->hasSubcube(*cube1));
    REQUIRE(cube0->hasSubcube(*cube2));
    REQUIRE(!(cube0->hasSubcube(*cube3)));
    REQUIRE(cube0->hasSubcube(*cube4));
    REQUIRE(cube0->hasSubcube(*cube5));
    REQUIRE(cube1->hasSubcube(*cube5));
    REQUIRE(cube2->hasSubcube(*cube4));
    REQUIRE(cube3->hasSubcube(*cube5));


    Cubeset cubeset{linkset, seqLens};
    std::cout << "Observed: " << cubeset.subcubeMap() << std::endl;
    std::cout << "Expected: " << expectedSubcubes << std::endl;
    REQUIRE(hopscotchEqual(cubeset.subcubeMap(), expectedSubcubes));
}
/*
TEST_CASE("Test Cube Scoring") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Cube Scoring" << std::endl;
    auto config = customConfiguration(CubeScoreParameter(5),
                                      CubeScoreThreshold(ULLONG_MAX),
                                      Genome1("genome0"),
                                      Genome2("genome1"),
                                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                         SpacedSeedMaskCollection::Span{5},
                                                                                                         SpacedSeedMaskCollection::SeedSetSize{1})},
                                      OldCubeScore(true),
                                      PerformGeometricHashing(true),
                                      TileSize(10));
    auto idMap = std::make_shared<IdentifierMapping>("genome0");
    idMap->queryGenomeID("genome1");
    idMap->querySequenceID("sequence0", "genome0");
    idMap->querySequenceID("sequence1", "genome1");

    SECTION("2D Case") {
        auto linkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMap);

        auto link0 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-0
                                                         KmerOccurrence(1, 1, 1, false, "ACGT")}, 4);
        linkset->addLink(link0);
        auto cube0 = std::make_shared<Cube>(*link0, config->tileSize()); // score 1, single Link

        auto link10_0 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-1
                                                            KmerOccurrence(1, 1, 11, false, "ACGT")}, 4);
        auto link10_1 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 10, false, "ACGT"),
                                                            KmerOccurrence(1, 1, 21, false, "ACGT")}, 4);
        linkset->addLink(link10_0);
        linkset->addLink(link10_1);
        auto cube1 = std::make_shared<Cube>(*link10_0, config->tileSize()); // score 2^2, two Links on same diagonal

        auto link20_0 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-2
                                                            KmerOccurrence(1, 1, 21, false, "ACGT")}, 4);
        auto link20_1 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),
                                                            KmerOccurrence(1, 1, 21, false, "ACGT")}, 4);
        linkset->addLink(link20_0);
        linkset->addLink(link20_1);
        auto cube2 = std::make_shared<Cube>(*link20_0, config->tileSize()); // score 1, single Link with count 2

        auto link30_0 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-3, bin 0
                                                            KmerOccurrence(1, 1, 31, false, "ACGT")}, 4);
        auto link30_1 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 1, false, "ACGT"),        // bin 1
                                                            KmerOccurrence(1, 1, 33, false, "ACGT")}, 4);
        auto link30_2 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 2, false, "ACGT"),
                                                            KmerOccurrence(1, 1, 34, false, "ACGT")}, 4);
        auto link30_3 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 3, false, "ACGT"),        // bin 2
                                                            KmerOccurrence(1, 1, 37, false, "ACGT")}, 4);
        linkset->addLink(link30_0);
        linkset->addLink(link30_1);
        linkset->addLink(link30_2);
        linkset->addLink(link30_3);
        auto cube3 = std::make_shared<Cube>(*link30_0, config->tileSize()); // score 6, one bin 0, two Links on same diagonal in bin 1, one bin 2

        auto link40_0 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-4, bin 0
                                                            KmerOccurrence(1, 1, 41, false, "ACGT")}, 4);
        auto link40_1 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 1, false, "ACGT"),        // bin 1
                                                            KmerOccurrence(1, 1, 43, false, "ACGT")}, 4);
        auto link40_2 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 2, false, "ACGT"),
                                                            KmerOccurrence(1, 1, 44, false, "ACGT")}, 4);
        auto link40_3 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 2, false, "ACGT"),
                                                            KmerOccurrence(1, 1, 44, false, "ACGT")}, 4);
        auto link40_4 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 3, false, "ACGT"),        // bin 2
                                                            KmerOccurrence(1, 1, 47, false, "ACGT")}, 4);
        auto link50_1 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 3, false, "ACGT"),        // --> now goes into cube 4 as diag +8, bin 4
                                                            KmerOccurrence(1, 1, 52, false, "ACGT")}, 4);
        linkset->addLink(link40_0);
        linkset->addLink(link40_1);
        linkset->addLink(link40_2);
        linkset->addLink(link40_3);
        linkset->addLink(link40_4);
        auto cube4 = std::make_shared<Cube>(*link40_0, config->tileSize()); // score 7, one bin 0, two Links on same diagonal, one of them with count 2 in bin 1, one bin 2, one in bin 4

        auto link50_0 = LinkPtr(std::vector<KmerOccurrence>{KmerOccurrence(0, 0, 0, false, "ACGT"),    // cube 0-5, diag -1
                                                            KmerOccurrence(1, 1, 51, false, "ACGT")}, 4);

        linkset->addLink(link50_0);
        linkset->addLink(link50_1);
        auto cube5 = std::make_shared<Cube>(*link50_0, config->tileSize()); // score 1
    }
}
*/
TEST_CASE("Test New Cube Scoring") {
    std::cout << "[INFO] -- [TEST CASE] -- Test New Cube Scoring" << std::endl;
    ConfigBuilder cbuilder;
    cbuilder.set(CubeLengthCutoff(600),
                 CubeScoreNormalizationParameter(3),
                 CubeScoreParameter(20), // nsubtiles
                 CubeScoreThreshold(0),
                 Hasse(true),
                 MatchLimit(SIZE_MAX),
                 OccurrencePerGenomeMax(SIZE_MAX),
                 OccurrencePerGenomeMin(1),
                 TileSize(100));
    auto config = cbuilder.makeConfig();
    auto idMap = std::make_shared<IdentifierMapping>("genome0");
    idMap->queryGenomeID("genome1");
    idMap->queryGenomeID("genome2");
    idMap->querySequenceID("sequence0", "genome0");
    idMap->querySequenceID("sequence1", "genome1");
    idMap->querySequenceID("sequence2", "genome2");
    auto seqLens = std::make_shared<tsl::hopscotch_map<size_t, size_t>>();
    (*seqLens)[0] = 1000;
    (*seqLens)[1] = 1500;
    (*seqLens)[2] = 700;
    size_t nPossibleLinks = 1; for (auto&& elem : *seqLens) { nPossibleLinks *= elem.second; }

    // make sure idMap is as expected
    REQUIRE(idMap->queryGenomeIDConst("genome0") == 0);
    REQUIRE(idMap->queryGenomeIDConst("genome1") == 1);
    REQUIRE(idMap->queryGenomeIDConst("genome2") == 2);
    REQUIRE(idMap->querySequenceName(0) == "sequence0");
    REQUIRE(idMap->querySequenceName(1) == "sequence1");
    REQUIRE(idMap->querySequenceName(2) == "sequence2");

    auto score = [&config](std::vector<double> const & bucketCounts, double p, double lambda, double n) {
        double score = 0;
        for (auto&& h : bucketCounts) {
            score += std::pow(h, p);
        }
        score = std::pow(score, 1./p);
        score /= (lambda * std::pow(n, 1./p));
        return score;
    };
    auto scoreEqual = [](double s1, double s2) {
        double f = 1.;
        while (s1*f < 1000.) { f *= 10; } // find factor to shift scores to at least 4 digits
        auto s1scale = std::floor(s1 * f)/f;
        auto s2scale = std::floor(s2 * f)/f;
        if (s1scale != s2scale) { std::cerr << "[DEBUG] -- scoreEqual -- " << s1 << " != " << s2 <<std::endl; }
        return s1scale == s2scale;
    };

    auto linkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMap);
    SECTION("2D Case") {
        // 20 subtiles, tilesize 100 -> i.e. diags [0,4], [5,9], [10,14], ..., [95, 99]
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,""},KmerOccurrence{1,1,0,false,""}}, 5); //diag 0, subtile 0
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,5,false,""},KmerOccurrence{1,1,5,false,""}}, 5);
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,10,false,""},KmerOccurrence{1,1,15,false,""}}, 5); // diag 5, subtile 1
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,15,false,""},KmerOccurrence{1,1,20,false,""}}, 5);
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,20,false,""},KmerOccurrence{1,1,118,false,""}}, 5); // diag 98, subtile 19

        // directly performs cube creation, cube scoring
        Cubeset cubeset{linkset, seqLens};

        REQUIRE(cubeset.cubeMap().size() == 1);
        auto cube = cubeset.cubeMap().begin();
        REQUIRE((cube->second).size() == 5);
        for (auto&& elem : cubeset.cubeMap()) {
            auto geo = elem.first->geoProperties(config->tileSize(), config->cubeScoreParameterChunks(), *seqLens);
            REQUIRE(geo == Cube::GeoProperties{1000, 0, 0});
        }
        REQUIRE(cubeset.computeCubeScore(cube->first) == score(std::vector<double>{2.,2.,1.},
                                                               config->cubeScoreNormalizationParameter(),
                                                               5.*1000.*5./(double)(nPossibleLinks), // chunkvol * nlinks / npossible
                                                               20.)); // b = 20
    }

    SECTION("3D Case") {
        // 20 subtiles, tilesize 100 -> 16 subtiles, 4 per side
        // i.e. diags [0,24], [25,49], [50,74], [75,99] (pairwise seq0 -- seq1/2)
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,""},
                                                         KmerOccurrence{1,1,0,false,""},        // diag 0
                                                         KmerOccurrence{2,2,0,false,""}}, 5);   // diag 0, subcube 0
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,5,false,""},
                                                         KmerOccurrence{1,1,55,false,""},       // diag 2
                                                         KmerOccurrence{2,2,5,false,""}}, 5);   // diag 0, subcube 2
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,5,false,""},
                                                         KmerOccurrence{1,1,55,false,""},       // diag 2
                                                         KmerOccurrence{2,2,55,false,""}}, 5);  // diag 2, subcube 10
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,10,false,""},
                                                         KmerOccurrence{1,1,35,false,""},       // diag 1
                                                         KmerOccurrence{2,2,60,false,""}}, 5);  // diag 2, subcube 9
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,15,false,""},
                                                         KmerOccurrence{1,1,40,false,""},       // diag 1
                                                         KmerOccurrence{2,2,65,false,""}}, 5);  // diag 2, subcube 9
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,20,false,""},
                                                         KmerOccurrence{1,1,20,false,""},       // diag 0
                                                         KmerOccurrence{2,2,45,false,""}}, 5);  // diag 1, subcube 4

        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,""},KmerOccurrence{1,1,0,false,""}}, 5); //diag 2, subtile 2
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,10,false,""},KmerOccurrence{2,2,35,false,""}}, 5); // diag 1, subtile 4

        auto linkfrac = 8./(double)(nPossibleLinks);
        // directly performs cube creation, cube scoring
        Cubeset cubeset{linkset, seqLens};

        std::shared_ptr<Cube const> cube;
        REQUIRE(cubeset.cubeMap().size() == 3);
        for (auto&& lcube : cubeset.cubeMap()) {
            if (lcube.first->dimensionality() == 3) {
                REQUIRE(lcube.second.size() == 6);
                cube = lcube.first;
                auto geo = cube->geoProperties(config->tileSize(), config->cubeScoreParameterChunks(), *seqLens);
                REQUIRE(geo == Cube::GeoProperties{700, 0, 0});
                REQUIRE(scoreEqual(cubeset.computeCubeScore(cube),
                                   score(std::vector<double>{1.,1.,1.,2.,1.},
                                         config->cubeScoreNormalizationParameter(),
                                         25.*25.*700.*linkfrac, // chunkvol * nlinks / npossible
                                         16.))); // 16 subtiles
            } else {
                REQUIRE(lcube.second.size() == 1);
                cube = lcube.first;
                if (cube->tiledistance(1).genome() == 1) {
                    auto geo = cube->geoProperties(config->tileSize(), config->cubeScoreParameterChunks(), *seqLens);
                    REQUIRE(geo == Cube::GeoProperties{1000, 0, 0});
                    REQUIRE(scoreEqual(cubeset.computeCubeScore(cube),
                                       score(std::vector<double>{1.},
                                             config->cubeScoreNormalizationParameter(),
                                             5.*1000.*linkfrac, // chunkvol * nlinks / npossible
                                             20.))); // b = 20
                } else {
                    auto geo = cube->geoProperties(config->tileSize(), config->cubeScoreParameterChunks(), *seqLens);
                    REQUIRE(geo == Cube::GeoProperties{700, 0, 0});
                    REQUIRE(scoreEqual(cubeset.computeCubeScore(cube),
                                       score(std::vector<double>{1.},
                                             config->cubeScoreNormalizationParameter(),
                                             5.*700.*linkfrac, // chunkvol * nlinks / npossible
                                             20.))); // b = 20
                }
            }
        }
    }

    SECTION("3D Case With Chunks") {
        cbuilder.set(CubeScoreParameterChunks(700));
        config = cbuilder.makeConfig();
        linkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMap);
        // 20 subtiles, tilesize 100 -> 16 subtiles, 4 per side
        // i.e. diags [0,24], [25,49], [50,74], [75,99] (pairwise seq0 -- seq1/2)
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,""},
                                                         KmerOccurrence{1,1,0,false,""},        // diag 0
                                                         KmerOccurrence{2,2,0,false,""}}, 5);   // diag 0, subcube 0, chunk 0
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,325,false,""},
                                                         KmerOccurrence{1,1,375,false,""},      // diag 2
                                                         KmerOccurrence{2,2,325,false,""}}, 5); // diag 0, subcube 2, chunk 1
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,350,false,""},
                                                         KmerOccurrence{1,1,400,false,""},      // diag 2
                                                         KmerOccurrence{2,2,350,false,""}}, 5); // diag 0, subcube 2, chunk 1
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,325,false,""},
                                                         KmerOccurrence{1,1,375,false,""},      // diag 2
                                                         KmerOccurrence{2,2,375,false,""}}, 5); // diag 2, subcube 10, chunk 1
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,10,false,""},
                                                         KmerOccurrence{1,1,35,false,""},       // diag 1
                                                         KmerOccurrence{2,2,60,false,""}}, 5);  // diag 2, subcube 9, chunk 0
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,615,false,""},
                                                         KmerOccurrence{1,1,640,false,""},      // diag 1
                                                         KmerOccurrence{2,2,665,false,""}}, 5); // diag 2, subcube 9, chunk 1
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,20,false,""},
                                                         KmerOccurrence{1,1,20,false,""},       // diag 0
                                                         KmerOccurrence{2,2,45,false,""}}, 5);  // diag 1, subcube 4, chunk 0

        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,1000,false,""},KmerOccurrence{1,1,1050,false,""}}, 5); // diag 2, subtile 2, chunk 3
        linkset->createLinks(std::vector<KmerOccurrence>{KmerOccurrence{0,0,10,false,""},KmerOccurrence{2,2,35,false,""}}, 5);     // diag 1, subtile 4, chunk 0
        auto linkfrac = 9./(double)(nPossibleLinks);

        // directly performs cube creation, cube scoring
        Cubeset cubeset{linkset, seqLens};
        std::shared_ptr<Cube const> cube;
        REQUIRE(cubeset.cubeMap().size() == 3);
        for (auto&& lcube : cubeset.cubeMap()) {
            if (lcube.first->dimensionality() == 3) {
                REQUIRE(lcube.second.size() == 7);
                cube = lcube.first;
                auto geo = cube->geoProperties(config->tileSize(), config->cubeScoreParameterChunks(), *seqLens);
                REQUIRE(geo == Cube::GeoProperties{700, 0, 3});
                REQUIRE(scoreEqual(cubeset.computeCubeScore(cube),
                                   score(std::vector<double>{1.,2.,1.,1.,1.,1.},
                                         config->cubeScoreNormalizationParameter(),
                                         25.*25.*700.*linkfrac, // chunkvol * nlinks / npossible
                                         64.))); // 16 subtiles * 4 chunks
            } else {
                REQUIRE(lcube.second.size() == 1);
                cube = lcube.first;
                if (cube->tiledistance(1).genome() == 1) {
                    auto geo = cube->geoProperties(config->tileSize(), config->cubeScoreParameterChunks(), *seqLens);
                    REQUIRE(geo == Cube::GeoProperties{1000, 0, 2});
                    REQUIRE(scoreEqual(cubeset.computeCubeScore(cube),
                                       score(std::vector<double>{1.},
                                             config->cubeScoreNormalizationParameter(),
                                             5.*700.*linkfrac, // chunkvol * nlinks / npossible
                                             60.))); // b=20 * 3 chunks
                } else {
                    auto geo = cube->geoProperties(config->tileSize(), config->cubeScoreParameterChunks(), *seqLens);
                    REQUIRE(geo == Cube::GeoProperties{700, 0, 1});
                    REQUIRE(scoreEqual(cubeset.computeCubeScore(cube),
                                       score(std::vector<double>{1.},
                                             config->cubeScoreNormalizationParameter(),
                                             5.*700.*linkfrac, // chunkvol * nlinks / npossible
                                             40.))); // b=20 * 2 chunks
                }
            }
        }
    }
}

TEST_CASE("Test Link Grouping") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Link Grouping" << std::endl;
    auto link1 = Link();
    link1.insertOccurrence(0, 0, 0, false, "ACGT");
    link1.insertOccurrence(1, 2, 0, false, "ACGT");
    link1.insertOccurrence(2, 4, 10, false, "ACGT");

    auto link2 = Link();  // same cube
    link2.insertOccurrence(0, 0, 10, false, "ACGT");
    link2.insertOccurrence(1, 2, 10, false, "ACGT");
    link2.insertOccurrence(2, 4, 20, false, "ACGT");

    auto link3 = Link();  // different cube
    link3.insertOccurrence(0, 0, 0, false, "ACGT");
    link3.insertOccurrence(1, 1, 0, false, "ACGT");
    link3.insertOccurrence(2, 4, 10, false, "ACGT");

    auto config = customConfiguration(CubeScoreParameter(3),
                                      CubeScoreThreshold(1),
                                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                         SpacedSeedMaskCollection::Span{5},
                                                                                                         SpacedSeedMaskCollection::SeedSetSize{1})},
                                      OldCubeScore(false),
                                      PerformGeometricHashing(true),
                                      TileSize(1));
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

    auto linkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMapping);
    linkset->addLink(LinkPtr(link1));
    linkset->addLink(LinkPtr(link1));
    linkset->addLink(LinkPtr(link2));
    linkset->addLink(LinkPtr(link3));
    REQUIRE(linkset->linkset().find(LinkPtr(link1)) != linkset->linkset().end());
    REQUIRE(linkset->linkset().find(LinkPtr(link2)) != linkset->linkset().end());
    REQUIRE(linkset->linkset().find(LinkPtr(link3)) != linkset->linkset().end());
    REQUIRE(linkset->linkset().at(LinkPtr(link1)) == 2); // Link count
    REQUIRE(linkset->linkset().at(LinkPtr(link2)) == 1);
    REQUIRE(linkset->linkset().at(LinkPtr(link3)) == 1);
}
/*
TEST_CASE("Test Link Count of Cubeset") {
    auto sameSequences = [](Cube const & cube, LinkPtr const & link) {
        if (cube.dimensionality() == link.dimensionality()) {
            for (size_t i = 0; i < cube.dimensionality(); ++i) {
                if (cube.tiledistance(i).genome() != link.genome(i)
                        || cube.tiledistance(i).sequence() != link.sequence(i)
                        || cube.tiledistance(i).reverse() != link.reverse(i)) {
                    return false;
                }
            }
            return true;
        } else {
            return false;
        }
    };

    std::cout << "[INFO] -- [TEST CASE] -- Test Link Count of Cubeset" << std::endl;
    auto data = LoadTestdataMetagraph();
    auto linksetFasta = data.runLinksetCreationFasta();
    auto expectedLinksFasta = data.expectedLinksFasta(*data.idMapGraph());
    REQUIRE(hopscotchEqual(linksetFasta->linkset(), expectedLinksFasta));
    Cubeset observedCubeset(linksetFasta, data.expectedSequenceLengthsExact());
    //observedCubeset.groupLinksInCubes();
    REQUIRE(hopscotchEqual(observedCubeset.cubeMap(), data.expectedCubesFasta(*data.idMapGraph())));
    tsl::hopscotch_map<std::shared_ptr<Cube const>, size_t, CubePtrHash, CubePtrEqual> expectedLinkcount;
    for (auto&& elem : observedCubeset.cubeMap()) { expectedLinkcount[elem.first] = 0; } // predefine cubes
    // count links by hand
    for (auto&& linkelem : expectedLinksFasta) {
        auto& link = linkelem.first;
        for (auto&& cubeelem : expectedLinkcount) {
            auto& cube = cubeelem.first;
            if (sameSequences(*cube, link)) {
                expectedLinkcount.at(cube) += 1;
            }
        }
    }
    // test link counts
    for (auto&& elem : expectedLinkcount) {
        REQUIRE(observedCubeset.linkCount(elem.first) == elem.second);
    }
}
*/
TEST_CASE("Test Expected Cubes with Metagraph") {
    auto roundScoreMap = [](Cubeset::ScoreMapType const & scoreToCube) {
        tsl::hopscotch_map<double, std::set<Cube>> roundedMap;
        std::vector<double> scores;
        auto minscore = scoreToCube.begin()->first;
        for (auto&& elem : scoreToCube) { minscore = std::min(minscore, elem.first); }
        double f = 1.;
        if (minscore < 1000.) {
            while (minscore*f < 1000.) { f *= 10; } // find factor to keep at least 4 digits
        } else {
            while (minscore*f > 1000.) { f /= 10; } // find factor to reduce to 4 digits
        }
        for (auto&& elem : scoreToCube) {
            auto s = std::floor(elem.first * f)/f;
            for (auto& cubeptr : elem.second) { roundedMap[s].insert(*cubeptr); }
            //roundedMap[s].insert(roundedMap[s].end(), elem.second.begin(), elem.second.end());
        }
        return roundedMap;
    };

    std::cout << "[INFO] -- [TEST CASE] -- Test Expected Cubes with Metagraph" << std::endl;
    auto data = LoadTestdataMetagraph();
    auto cubesetFasta = data.runCubesetCreationFasta();
    auto cubesetGraph = data.runCubesetCreationGraph();
    //cubesetFasta.groupLinksInCubes();
    //cubesetGraph.groupLinksInCubes();
    auto expectedCubesFasta = data.expectedCubesFasta(*data.idMapGraph());
    auto expectedCubesGraph = data.expectedCubesGraph(*data.idMapGraph());
    auto expectedCubeScoresFasta = data.expectedCubeScoresFasta(*data.idMapGraph());
    auto expectedCubeScoresGraph = data.expectedCubeScoresGraph(*data.idMapGraph());
    REQUIRE(hopscotchEqual(expectedCubesFasta, cubesetFasta.cubeMap()));
    REQUIRE(hopscotchEqual(expectedCubesGraph, cubesetGraph.cubeMap()));
std::cout << "[DEBUG] -- expectedScores\n" << roundScoreMap(expectedCubeScoresFasta) << std::endl << std::endl << std::endl;
std::cout << "[DEBUG] -- scoreToCube\n" << roundScoreMap(cubesetFasta.scoreToCube()) << std::endl << std::endl << std::endl;
    REQUIRE(hopscotchEqual(roundScoreMap(expectedCubeScoresFasta), roundScoreMap(cubesetFasta.scoreToCube())));
    REQUIRE(hopscotchEqual(roundScoreMap(expectedCubeScoresGraph), roundScoreMap(cubesetGraph.scoreToCube())));
    // subcube link inclusion
    auto expectedCubesWithSubcubeLinksFasta = data.expectedCubesWithSubcubeLinksFasta(*data.idMapGraph());
    auto expectedCubesWithSubcubeLinksGraph = data.expectedCubesWithSubcubeLinksGraph(*data.idMapGraph());
    tsl::hopscotch_map<std::shared_ptr<Cube const>,
                       tsl::hopscotch_set<LinkPtr, LinkPtrHash>,
                       CubePtrHash, CubePtrEqual> observedCubesWithSubcubeLinksFasta;
    tsl::hopscotch_map<std::shared_ptr<Cube const>,
                       tsl::hopscotch_set<LinkPtr, LinkPtrHash>,
                       CubePtrHash, CubePtrEqual> observedCubesWithSubcubeLinksGraph;
    for (auto&& cube : cubesetFasta.cubeMap()) {
        observedCubesWithSubcubeLinksFasta.insert({cube.first, cubesetFasta.linksIncludingSubcubes(cube.first)});
        //REQUIRE(expectedCubesWithSubcubeLinksFasta.at(cube.first) == cubesetFasta.linksIncludingSubcubes(cube.first));
    }
    for (auto&& cube : cubesetGraph.cubeMap()) {
        observedCubesWithSubcubeLinksGraph.insert({cube.first, cubesetGraph.linksIncludingSubcubes(cube.first)});
        //REQUIRE(expectedCubesWithSubcubeLinksGraph.at(cube.first) == cubesetGraph.linksIncludingSubcubes(cube.first));
    }
    REQUIRE(hopscotchEqual(expectedCubesWithSubcubeLinksFasta, observedCubesWithSubcubeLinksFasta));
    REQUIRE(hopscotchEqual(expectedCubesWithSubcubeLinksGraph, observedCubesWithSubcubeLinksGraph));
}

TEST_CASE("Test Region Tuple Extraction") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Region Tuple Extraction" << std::endl;
    auto data = LoadTestdataMetagraph();
    auto tuplesFasta = data.runRegionTupleExtractionFasta();
    auto tuplesGraph = data.runRegionTupleExtractionGraph();
    auto expectedTuplesFasta = data.expectedRegionTuplesFasta(*data.idMapGraph());
    auto expectedTuplesGraph = data.expectedRegionTuplesGraph(*data.idMapGraph());
    REQUIRE(hopscotchEqual(expectedTuplesFasta, tuplesFasta));
    REQUIRE(hopscotchEqual(expectedTuplesGraph, tuplesGraph));
//std::cout << "[DEBUG] -- testCubeRelated -- fasta tuples: " << tuplesFasta << std::endl;
//std::cout << "[DEBUG] -- testCubeRelated -- graph tuples: " << tuplesGraph << std::endl;
    for (auto& elem : tuplesFasta) {
        size_t prevStart = SIZE_MAX;
        size_t prevEnd = SIZE_MAX;
        for (auto& tuple : elem.second) { // tuples in descending order
            auto tstart = tuple.at(0);
            auto tend = tuple.at(1);
            REQUIRE(tstart < prevStart);
            REQUIRE(tend < prevEnd); // strictly sorted
            REQUIRE(tend < prevStart); // no overlap
            prevStart = tstart;
            prevEnd = tend;
        }
    }
    for (auto& elem : tuplesGraph) {
        size_t prevStart = SIZE_MAX;
        size_t prevEnd = SIZE_MAX;
        for (auto& tuple : elem.second) { // tuples in descending order
            auto tstart = tuple.at(0);
            auto tend = tuple.at(1);
            REQUIRE(tstart < prevStart);
            REQUIRE(tend < prevEnd); // strictly sorted
            REQUIRE(tend < prevStart); // no overlap
            prevStart = tstart;
            prevEnd = tend;
        }
    }
}

TEST_CASE("Test PrefilterCubeset") {
    std::cout << "[INFO] -- [TEST CASE] -- Test PrefilterCubeset" << std::endl;
    auto data = LoadTestdataMetagraph();
    auto seedMap = data.runSeedExtractionFasta();
    auto idMap = data.idMapFasta();
    auto seqlens = data.expectedSequenceLengthsExact();
//std::cout << std::endl << std::endl << "[DEBUG] -- 000 -- seedMap " << seedMap->seedMap() << std::endl << std::endl << std::endl;
    ConfigBuilder configBuilder{};
    configBuilder.set(CubeScoreThreshold{0},
                      Genome1{data.configFasta()->genome1()},
                      Genome2{data.configFasta()->genome2()},
                      Hasse(true),
                      MaskCollectionPtr{data.configFasta()->maskCollection()},
                      NThreads{3},
                      MatchLimit{ULLONG_MAX},
                      OptimalSeed{true},
                      PerformGeometricHashing{true},
                      PreHasse{true},
                      PreLinkThreshold{0},
                      PreMaskCollectionPtr{data.configFasta()->maskCollection()},
                      PreOptimalSeed{true},
                      TileSize{500});

    SECTION("No Link Threshold") {
//std::cout << std::endl << std::endl << "[DEBUG] -- 001" << std::endl << std::endl << std::endl;
        auto config = configBuilder.makeConfig();
        auto linkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMap);
        linkset->createLinks(*seedMap);
        Cubeset cubeSet(linkset, seqlens);

        PrefilterCubeset prefilterCubeset(*seedMap, idMap, seqlens, config);
        REQUIRE(cubeSet.cubeMap().size() == prefilterCubeset.relevantCubeSet().size());
        REQUIRE(cubeSet.cubeMap().size() == prefilterCubeset.cubeMap().size());
        for (auto&& elem : cubeSet.cubeMap()) {
            REQUIRE(prefilterCubeset.relevantCubeSet().find(elem.first) != prefilterCubeset.relevantCubeSet().end());
            REQUIRE(prefilterCubeset.cubeMap().find(elem.first) != prefilterCubeset.cubeMap().end());
            size_t count = 0;
            for (auto&& link : elem.second) {
                count += cubeSet.underlyingLinkset()->linkCount(link);
            }
            REQUIRE(count == prefilterCubeset.cubeMap().at(elem.first));
        }
        for (auto&& elem : prefilterCubeset.relevantCubeSet()) {
            REQUIRE(cubeSet.cubeMap().find(elem) != cubeSet.cubeMap().end());
        }
        for (auto&& elem : prefilterCubeset.cubeMap()) {
            REQUIRE(cubeSet.cubeMap().find(elem.first) != cubeSet.cubeMap().end());
        }
    }

    SECTION("With Link Threshold") {
//std::cout << std::endl << std::endl << "[DEBUG] -- 002" << std::endl << std::endl << std::endl;
        size_t threshold = 3;
        configBuilder.set(PreLinkThreshold{threshold});
        auto config = configBuilder.makeConfig();
        auto linkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMap);
        linkset->createLinks(*seedMap);
        Cubeset cubeSet(linkset, seqlens);

        PrefilterCubeset prefilterCubeset(*seedMap, idMap, seqlens, config);
        REQUIRE(cubeSet.cubeMap().size() > prefilterCubeset.relevantCubeSet().size());
        REQUIRE(cubeSet.cubeMap().size() == prefilterCubeset.cubeMap().size());
        for (auto&& elem : cubeSet.cubeMap()) {
            REQUIRE(prefilterCubeset.cubeMap().find(elem.first) != prefilterCubeset.cubeMap().end());
            size_t count = 0;
            for (auto&& link : elem.second) {
                count += cubeSet.underlyingLinkset()->linkCount(link);
            }
            REQUIRE(count == prefilterCubeset.cubeMap().at(elem.first));
            if (count >= threshold) {
                REQUIRE(prefilterCubeset.relevantCubeSet().find(elem.first) != prefilterCubeset.relevantCubeSet().end());
            }
        }
        for (auto&& elem : prefilterCubeset.relevantCubeSet()) {
            REQUIRE(cubeSet.cubeMap().find(elem) != cubeSet.cubeMap().end());
        }
        for (auto&& elem : prefilterCubeset.cubeMap()) {
            REQUIRE(cubeSet.cubeMap().find(elem.first) != cubeSet.cubeMap().end());
        }
    }
}



auto pGHExpectedSeeds(IdentifierMapping const & idMap, fs::path file) {
    auto expectedSeeds = typename SeedMap<TwoBitKmerDataShort>::SeedMapType{};
    // load json
    std::ifstream is(file);
    nlohmann::json j;
    is >> j;
    // fill expectedSeeds
    for (auto&& elem : j.items()) {
        auto seed = TwoBitKmer<TwoBitKmerDataShort>(elem.key());
        std::vector<KmerOccurrence> occurrences{};
        for (auto&& occ : elem.value()) {
            if (occ[0].get<std::string>() == "hg38") {
                occurrences.emplace_back(0, idMap.querySequenceIDConst(occ[1].get<std::string>(), "hg38"), occ[2].get<size_t>(), false, elem.key());
            } else {
                occurrences.emplace_back(1, idMap.querySequenceIDConst(occ[1].get<std::string>(), "mm10"), occ[2].get<size_t>(), false, elem.key());
            }
            std::sort(occurrences.begin(), occurrences.end());
        }
        expectedSeeds[seed].emplace_back(occurrences);
    }

    return expectedSeeds;
}

TEST_CASE("Test PrefilterCubeset on Real Sequences") {
    std::cout << "[INFO] -- [TEST CASE] -- Test PrefilterCubeset on Real Sequences" << std::endl;
    auto filepath = testfilesBasepath() / "pGH";
    ConfigBuilder configBuilder{};
    configBuilder.set(CubeScoreThreshold{0},
                      Genome1{"hg38"},
                      Genome2{"mm10"},
                      Hasse(true),
                      InputFiles{std::vector<std::string>{(filepath/"hg38.fa").string(), (filepath/"mm10.fa").string()}},
                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{12})},
                      NThreads{3},
                      MatchLimit{ULLONG_MAX},
                      OccurrencePerSequenceMax{ULLONG_MAX},
                      OptimalSeed{true},
                      PerformGeometricHashing{true},
                      PreHasse{true},
                      PreLinkThreshold{0},
                      PreMaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{12})},
                      PreOptimalSeed{true},
                      Redmask{false},
                      TileSize{10000},
                      Thinning{1});
    auto config = configBuilder.makeConfig();

    auto fastaCollection = std::make_shared<FastaCollection const>(config);
    auto idMap = std::make_shared<IdentifierMapping>("hg38");
    auto seqLens = std::make_shared<tsl::hopscotch_map<size_t, size_t>>();
    fastaCollection->populateIdentifierMappingFromFastaCollection(*idMap);
    fastaCollection->fillSequenceLengths(*seqLens, *idMap);

    // assert ids and get sequence ids
    REQUIRE(idMap->numGenomes() == 2);
    REQUIRE(idMap->numSequences() == 14);
    REQUIRE(idMap->queryGenomeName(0) == "hg38");
    REQUIRE(idMap->queryGenomeName(1) == "mm10");

    SeedMap<TwoBitKmerDataShort> seedMap(config, idMap, PrefilterSeedMapTag{});
    seedMap.extractSeeds(fastaCollection);
    auto expectedSeedMap = pGHExpectedSeeds(*idMap, filepath/"pGhTestSeeds.json");
    auto& observedSeedMap = seedMap.seedMap();

    auto compareSeedMaps = [&]() {
        bool equal = true;
        if (observedSeedMap.size() != expectedSeedMap.size()) {
            std::cout << "[ERR] -- " << observedSeedMap.size() << " (obs) != " << expectedSeedMap.size() << " (exp)" << std::endl;
            equal = false;
        }
        for (auto&& elem : expectedSeedMap) {
            REQUIRE(observedSeedMap.find(elem.first) != observedSeedMap.end());
        }
        for (auto&& elem : observedSeedMap) {
            REQUIRE(expectedSeedMap.find(elem.first) != expectedSeedMap.end());
        }
        for (auto&& elem : expectedSeedMap) {
            auto& seed = elem.first;
            auto& occv = elem.second;
            if (observedSeedMap.find(seed) == observedSeedMap.end()) {
                std::cout << "[ERR] -- " << seed << " not in observed seed map" << std::endl;
                equal = false;
            } else {
                REQUIRE(occv.size() == 1);
                REQUIRE(observedSeedMap.at(seed).size() == 1);
                auto expOccs = occv.at(0);
                auto obsOccs = observedSeedMap.at(seed).at(0);
                std::sort(obsOccs.begin(), obsOccs.end());
                if (expOccs.size() != obsOccs.size()) {
                    std::cout << "[ERR] -- " << seed << " -- " << expOccs << " (exp) != " << obsOccs << " (obs)" << std::endl;
                    equal = false;
                } else {
                    auto veq = true;
                    for (size_t i = 0; i < expOccs.size(); ++i) {
                        if (!KmerOccurrence::equalSpot(expOccs.at(i), obsOccs.at(i))) { veq = false; }
                    }
                    if (!veq) {
                        std::cout << "[ERR] -- " << seed << " -- " << expOccs << " (exp) != " << obsOccs << " (obs)" << std::endl;
                        equal = false;
                    }
                }
            }
        }
        return equal;
    };

    REQUIRE(compareSeedMaps());

    // reference cubeset from python
    tsl::hopscotch_map<std::shared_ptr<Cube const>, size_t,
                       CubePtrHash, CubePtrEqual> requiredRelevantCubeMap;

    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:33650593|genestart:33657087|rgenestart:6495|seqlen:29156|genelen:17377|chrlen:64444167|chr:chr20|strand:-|gid:ENSG00000125967|tid:ENST00000375238", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:154536355|genestart:154544399|rgenestart:8045|seqlen:32015|genelen:14492|chrlen:182113224|chr:chr2|strand:-|gid:ENSMUSG00000027489|tid:ENSMUSG00000027489", "mm10"),0,false} } )] = 448;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:33650593|genestart:33657087|rgenestart:6495|seqlen:29156|genelen:17377|chrlen:64444167|chr:chr20|strand:-|gid:ENSG00000125967|tid:ENST00000375238", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:154536355|genestart:154544399|rgenestart:8045|seqlen:32015|genelen:14492|chrlen:182113224|chr:chr2|strand:-|gid:ENSMUSG00000027489|tid:ENSMUSG00000027489", "mm10"),2,false} } )] = 9;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:33650593|genestart:33657087|rgenestart:6495|seqlen:29156|genelen:17377|chrlen:64444167|chr:chr20|strand:-|gid:ENSG00000125967|tid:ENST00000375238", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:154536355|genestart:154544399|rgenestart:8045|seqlen:32015|genelen:14492|chrlen:182113224|chr:chr2|strand:-|gid:ENSMUSG00000027489|tid:ENSMUSG00000027489", "mm10"),1,false} } )] = 8;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:33650593|genestart:33657087|rgenestart:6495|seqlen:29156|genelen:17377|chrlen:64444167|chr:chr20|strand:-|gid:ENSG00000125967|tid:ENST00000375238", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:154536355|genestart:154544399|rgenestart:8045|seqlen:32015|genelen:14492|chrlen:182113224|chr:chr2|strand:-|gid:ENSMUSG00000027489|tid:ENSMUSG00000027489", "mm10"),-1,false} } )] = 17;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:33650593|genestart:33657087|rgenestart:6495|seqlen:29156|genelen:17377|chrlen:64444167|chr:chr20|strand:-|gid:ENSG00000125967|tid:ENST00000375238", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:154536355|genestart:154544399|rgenestart:8045|seqlen:32015|genelen:14492|chrlen:182113224|chr:chr2|strand:-|gid:ENSMUSG00000027489|tid:ENSMUSG00000027489", "mm10"),-2,false} } )] = 11;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:33650593|genestart:33657087|rgenestart:6495|seqlen:29156|genelen:17377|chrlen:64444167|chr:chr20|strand:-|gid:ENSG00000125967|tid:ENST00000375238", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:154536355|genestart:154544399|rgenestart:8045|seqlen:32015|genelen:14492|chrlen:182113224|chr:chr2|strand:-|gid:ENSMUSG00000027489|tid:ENSMUSG00000027489", "mm10"),-3,false} } )] = 3;

    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:100863632|genestart:100872372|rgenestart:8741|seqlen:37561|genelen:23627|chrlen:248956422|chr:chr1|strand:-|gid:ENSG00000162694|tid:ENST00000370113", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:116001804|genestart:116007462|rgenestart:5659|seqlen:35117|genelen:21556|chrlen:160039680|chr:chr3|strand:+|gid:ENSMUSG00000027963|tid:ENSMUSG00000027963", "mm10"),-1,false} } )] = 172;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:100863632|genestart:100872372|rgenestart:8741|seqlen:37561|genelen:23627|chrlen:248956422|chr:chr1|strand:-|gid:ENSG00000162694|tid:ENST00000370113", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:116001804|genestart:116007462|rgenestart:5659|seqlen:35117|genelen:21556|chrlen:160039680|chr:chr3|strand:+|gid:ENSMUSG00000027963|tid:ENSMUSG00000027963", "mm10"),2,false} } )] = 4;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:100863632|genestart:100872372|rgenestart:8741|seqlen:37561|genelen:23627|chrlen:248956422|chr:chr1|strand:-|gid:ENSG00000162694|tid:ENST00000370113", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:116001804|genestart:116007462|rgenestart:5659|seqlen:35117|genelen:21556|chrlen:160039680|chr:chr3|strand:+|gid:ENSMUSG00000027963|tid:ENSMUSG00000027963", "mm10"),1,false} } )] = 12;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:100863632|genestart:100872372|rgenestart:8741|seqlen:37561|genelen:23627|chrlen:248956422|chr:chr1|strand:-|gid:ENSG00000162694|tid:ENST00000370113", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:116001804|genestart:116007462|rgenestart:5659|seqlen:35117|genelen:21556|chrlen:160039680|chr:chr3|strand:+|gid:ENSMUSG00000027963|tid:ENSMUSG00000027963", "mm10"),0,false} } )] = 217;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:100863632|genestart:100872372|rgenestart:8741|seqlen:37561|genelen:23627|chrlen:248956422|chr:chr1|strand:-|gid:ENSG00000162694|tid:ENST00000370113", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:116001804|genestart:116007462|rgenestart:5659|seqlen:35117|genelen:21556|chrlen:160039680|chr:chr3|strand:+|gid:ENSMUSG00000027963|tid:ENSMUSG00000027963", "mm10"),-2,false} } )] = 11;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:100863632|genestart:100872372|rgenestart:8741|seqlen:37561|genelen:23627|chrlen:248956422|chr:chr1|strand:-|gid:ENSG00000162694|tid:ENST00000370113", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:116001804|genestart:116007462|rgenestart:5659|seqlen:35117|genelen:21556|chrlen:160039680|chr:chr3|strand:+|gid:ENSMUSG00000027963|tid:ENSMUSG00000027963", "mm10"),-3,false} } )] = 10;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:100863632|genestart:100872372|rgenestart:8741|seqlen:37561|genelen:23627|chrlen:248956422|chr:chr1|strand:-|gid:ENSG00000162694|tid:ENST00000370113", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:116001804|genestart:116007462|rgenestart:5659|seqlen:35117|genelen:21556|chrlen:160039680|chr:chr3|strand:+|gid:ENSMUSG00000027963|tid:ENSMUSG00000027963", "mm10"),-4,false} } )] = 2;

    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),8,false} } )] = 10;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),6,false} } )] = 18;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),0,false} } )] = 37;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),1,false} } )] = 24;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),2,false} } )] = 23;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),12,false} } )] = 1;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-1,false} } )] = 216;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),11,false} } )] = 5;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),4,false} } )] = 15;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),5,false} } )] = 19;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),9,false} } )] = 7;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),7,false} } )] = 15;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),3,false} } )] = 26;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-10,false} } )] = 17;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),10,false} } )] = 4;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-4,false} } )] = 52;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-3,false} } )] = 43;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-2,false} } )] = 121;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-7,false} } )] = 24;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-11,false} } )] = 20;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-5,false} } )] = 27;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-6,false} } )] = 30;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-8,false} } )] = 9;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-9,false} } )] = 15;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-13,false} } )] = 19;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-12,false} } )] = 15;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-15,false} } )] = 2;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855", "mm10"),-14,false} } )] = 7;

    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:35607982|genestart:35615892|rgenestart:7911|seqlen:21825|genelen:5158|chrlen:64444167|chr:chr20|strand:+|gid:ENSG00000061656|tid:ENST00000374273", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:156058586|genestart:156065175|rgenestart:6590|seqlen:19939|genelen:5443|chrlen:182113224|chr:chr2|strand:+|gid:ENSMUSG00000038180|tid:ENSMUSG00000038180", "mm10"),0,false} } )] = 5;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:35607982|genestart:35615892|rgenestart:7911|seqlen:21825|genelen:5158|chrlen:64444167|chr:chr20|strand:+|gid:ENSG00000061656|tid:ENST00000374273", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:156058586|genestart:156065175|rgenestart:6590|seqlen:19939|genelen:5443|chrlen:182113224|chr:chr2|strand:+|gid:ENSMUSG00000038180|tid:ENSMUSG00000038180", "mm10"),-1,false} } )] = 312;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:35607982|genestart:35615892|rgenestart:7911|seqlen:21825|genelen:5158|chrlen:64444167|chr:chr20|strand:+|gid:ENSG00000061656|tid:ENST00000374273", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:156058586|genestart:156065175|rgenestart:6590|seqlen:19939|genelen:5443|chrlen:182113224|chr:chr2|strand:+|gid:ENSMUSG00000038180|tid:ENSMUSG00000038180", "mm10"),-2,false} } )] = 1;

    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),0,false} } )] = 12;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),3,false} } )] = 1;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),2,false} } )] = 4;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),1,false} } )] = 7;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),-1,false} } )] = 60;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),-2,false} } )] = 84;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),-4,false} } )] = 14;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),-3,false} } )] = 49;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),-5,false} } )] = 6;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824", "mm10"),-6,false} } )] = 2;

    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),1,false} } )] = 40;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),5,false} } )] = 149;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),10,false} } )] = 12;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),-1,false} } )] = 33;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),12,false} } )] = 5;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),9,false} } )] = 25;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),0,false} } )] = 32;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),4,false} } )] = 40;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),3,false} } )] = 35;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),6,false} } )] = 356;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),11,false} } )] = 6;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),7,false} } )] = 39;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),8,false} } )] = 33;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),2,false} } )] = 27;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),-2,false} } )] = 21;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),-3,false} } )] = 23;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),-4,false} } )] = 17;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),-5,false} } )] = 11;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447", "mm10"),-6,false} } )] = 8;

    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:123374587|genestart:123383773|rgenestart:9187|seqlen:43706|genelen:25586|chrlen:133275309|chr:chr12|strand:+|gid:ENSG00000183955|tid:ENST00000402868", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:124432181|genestart:124439930|rgenestart:7750|seqlen:36339|genelen:22379|chrlen:151834684|chr:chr5|strand:+|gid:ENSMUSG00000049327|tid:ENSMUSG00000049327", "mm10"),1,false} } )] = 7;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:123374587|genestart:123383773|rgenestart:9187|seqlen:43706|genelen:25586|chrlen:133275309|chr:chr12|strand:+|gid:ENSG00000183955|tid:ENST00000402868", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:124432181|genestart:124439930|rgenestart:7750|seqlen:36339|genelen:22379|chrlen:151834684|chr:chr5|strand:+|gid:ENSMUSG00000049327|tid:ENSMUSG00000049327", "mm10"),0,false} } )] = 15;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:123374587|genestart:123383773|rgenestart:9187|seqlen:43706|genelen:25586|chrlen:133275309|chr:chr12|strand:+|gid:ENSG00000183955|tid:ENST00000402868", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:124432181|genestart:124439930|rgenestart:7750|seqlen:36339|genelen:22379|chrlen:151834684|chr:chr5|strand:+|gid:ENSMUSG00000049327|tid:ENSMUSG00000049327", "mm10"),2,false} } )] = 3;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:123374587|genestart:123383773|rgenestart:9187|seqlen:43706|genelen:25586|chrlen:133275309|chr:chr12|strand:+|gid:ENSG00000183955|tid:ENST00000402868", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:124432181|genestart:124439930|rgenestart:7750|seqlen:36339|genelen:22379|chrlen:151834684|chr:chr5|strand:+|gid:ENSMUSG00000049327|tid:ENSMUSG00000049327", "mm10"),-1,false} } )] = 544;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:123374587|genestart:123383773|rgenestart:9187|seqlen:43706|genelen:25586|chrlen:133275309|chr:chr12|strand:+|gid:ENSG00000183955|tid:ENST00000402868", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:124432181|genestart:124439930|rgenestart:7750|seqlen:36339|genelen:22379|chrlen:151834684|chr:chr5|strand:+|gid:ENSMUSG00000049327|tid:ENSMUSG00000049327", "mm10"),-2,false} } )] = 6;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:123374587|genestart:123383773|rgenestart:9187|seqlen:43706|genelen:25586|chrlen:133275309|chr:chr12|strand:+|gid:ENSG00000183955|tid:ENST00000402868", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:124432181|genestart:124439930|rgenestart:7750|seqlen:36339|genelen:22379|chrlen:151834684|chr:chr5|strand:+|gid:ENSMUSG00000049327|tid:ENSMUSG00000049327", "mm10"),-3,false} } )] = 3;
    requiredRelevantCubeMap[std::make_shared<Cube const>( std::vector<Tiledistance>{ Tiledistance{0,(uint32_t)idMap->querySequenceIDConst("genome:hg38|seqstart:123374587|genestart:123383773|rgenestart:9187|seqlen:43706|genelen:25586|chrlen:133275309|chr:chr12|strand:+|gid:ENSG00000183955|tid:ENST00000402868", "hg38"),0,false}, Tiledistance{1,(uint32_t)idMap->querySequenceIDConst("genome:mm10|seqstart:124432181|genestart:124439930|rgenestart:7750|seqlen:36339|genelen:22379|chrlen:151834684|chr:chr5|strand:+|gid:ENSMUSG00000049327|tid:ENSMUSG00000049327", "mm10"),-4,false} } )] = 1;

    // actual cubeset
    PrefilterCubeset relevantCubeset{seedMap, idMap, seqLens, config};
    Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset(config, idMap);
    linkset.createLinks(seedMap);

    // re-code cube counting
    tsl::hopscotch_map<Link, size_t, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> recodeLinkset;
    tsl::hopscotch_map<std::shared_ptr<Cube const>, size_t,
                       CubePtrHash, CubePtrEqual> recodeRelevantCubeMap;
    //std::ofstream seedOut(filepath/"GHseeds.json");
    //seedOut << "[";
    //REQUIRE(seedOut.good());
    for (auto&& elem : observedSeedMap) {
        std::vector<KmerOccurrence> hgOccs;
        std::vector<KmerOccurrence> mmOccs;
        REQUIRE(elem.second.size() == 1);
        for (auto& occ : elem.second.at(0)) {
            REQUIRE(occ.genome() >= 0);
            REQUIRE(occ.genome() < 2);
            if (occ.genome() == 0) {
                hgOccs.emplace_back(occ);
            } else {
                mmOccs.emplace_back(occ);
            }
        }
        REQUIRE((hgOccs.size() + mmOccs.size()) == elem.second.at(0).size());
        if (hgOccs.size() > 0 && mmOccs.size() > 0) {
            for (auto&& hgOcc : hgOccs) {
                for (auto&& mmOcc : mmOccs) {
                    Link link(std::vector<KmerOccurrence>{hgOcc, mmOcc}, config->preMaskCollection()->span(0));
                    if (recodeLinkset.find(link) == recodeLinkset.end()) {
                        recodeLinkset[link] = 0;
                    }
                    recodeLinkset.at(link) += 1;
                    auto cube = std::make_shared<Cube>(link, config->tileSize());
                    if (recodeRelevantCubeMap.find(cube) == recodeRelevantCubeMap.end()) {
                        recodeRelevantCubeMap[cube] = 0;
                    }
                    recodeRelevantCubeMap.at(cube) += 1;

                    REQUIRE(cube->dimensionality() == 2);
                    auto delta = static_cast<long long>(mmOcc.position()) - static_cast<long long>(hgOcc.position());
                    auto tile = static_cast<long long>(std::floor(static_cast<double>(delta)/static_cast<double>(config->tileSize())));
                    REQUIRE(cube->tiledistance(1).distance() == tile);

                    //seedOut << "[[\"" << idMap->querySequenceName(hgOcc.sequence()) << "\", " << hgOcc.position() << "],";
                    //seedOut <<  "[\"" << idMap->querySequenceName(mmOcc.sequence()) << "\", " << mmOcc.position() << "]],\n";
                }
            }
        }
    }
    //seedOut << "]";

    REQUIRE(hopscotchEqual(recodeLinkset, linkset.linkset()));
    REQUIRE(hopscotchEqual(recodeRelevantCubeMap, relevantCubeset.cubeMap()));
    // check if required cubes are in unfiltered data
    for (auto&& elem : requiredRelevantCubeMap) {
        REQUIRE(recodeRelevantCubeMap.find(elem.first) != recodeRelevantCubeMap.end());
        REQUIRE(recodeRelevantCubeMap.at(elem.first) == elem.second);
    }

    // filter relevantCubeset.cubeMap() for only orthologous cubes
    tsl::hopscotch_map<std::string, std::string> orthology;
    orthology["genome:hg38|seqstart:33650593|genestart:33657087|rgenestart:6495|seqlen:29156|genelen:17377|chrlen:64444167|chr:chr20|strand:-|gid:ENSG00000125967|tid:ENST00000375238"] = "genome:mm10|seqstart:154536355|genestart:154544399|rgenestart:8045|seqlen:32015|genelen:14492|chrlen:182113224|chr:chr2|strand:-|gid:ENSMUSG00000027489|tid:ENSMUSG00000027489";
    orthology["genome:hg38|seqstart:100863632|genestart:100872372|rgenestart:8741|seqlen:37561|genelen:23627|chrlen:248956422|chr:chr1|strand:-|gid:ENSG00000162694|tid:ENST00000370113"] = "genome:mm10|seqstart:116001804|genestart:116007462|rgenestart:5659|seqlen:35117|genelen:21556|chrlen:160039680|chr:chr3|strand:+|gid:ENSMUSG00000027963|tid:ENSMUSG00000027963";
    orthology["genome:hg38|seqstart:114847033|genestart:114854803|rgenestart:7771|seqlen:156023|genelen:140568|chrlen:248956422|chr:chr1|strand:+|gid:ENSG00000198765|tid:ENST00000369518"] = "genome:mm10|seqstart:102810321|genestart:102818499|rgenestart:8179|seqlen:132474|genelen:117602|chrlen:160039680|chr:chr3|strand:-|gid:ENSMUSG00000027855|tid:ENSMUSG00000027855";
    orthology["genome:hg38|seqstart:35607982|genestart:35615892|rgenestart:7911|seqlen:21825|genelen:5158|chrlen:64444167|chr:chr20|strand:+|gid:ENSG00000061656|tid:ENST00000374273"] = "genome:mm10|seqstart:156058586|genestart:156065175|rgenestart:6590|seqlen:19939|genelen:5443|chrlen:182113224|chr:chr2|strand:+|gid:ENSMUSG00000038180|tid:ENSMUSG00000038180";
    orthology["genome:hg38|seqstart:66502723|genestart:66508003|rgenestart:5281|seqlen:56870|genelen:44542|chrlen:90338345|chr:chr16|strand:-|gid:ENSG00000166548|tid:ENST00000299697"] = "genome:mm10|seqstart:104217325|genestart:104226685|rgenestart:9361|seqlen:36971|genelen:21874|chrlen:129401213|chr:chr8|strand:-|gid:ENSMUSG00000035824|tid:ENSMUSG00000035824";
    orthology["genome:hg38|seqstart:108035287|genestart:108043094|rgenestart:7808|seqlen:64705|genelen:48769|chrlen:198295559|chr:chr3|strand:-|gid:ENSG00000196776|tid:ENST00000361309"] = "genome:mm10|seqstart:49793395|genestart:49800533|rgenestart:7139|seqlen:129559|genelen:114478|chrlen:98207768|chr:chr16|strand:+|gid:ENSMUSG00000055447|tid:ENSMUSG00000055447";
    orthology["genome:hg38|seqstart:123374587|genestart:123383773|rgenestart:9187|seqlen:43706|genelen:25586|chrlen:133275309|chr:chr12|strand:+|gid:ENSG00000183955|tid:ENST00000402868"] = "genome:mm10|seqstart:124432181|genestart:124439930|rgenestart:7750|seqlen:36339|genelen:22379|chrlen:151834684|chr:chr5|strand:+|gid:ENSMUSG00000049327|tid:ENSMUSG00000049327";
    tsl::hopscotch_map<std::shared_ptr<Cube const>, size_t,
                       CubePtrHash, CubePtrEqual> orthologRelevantCubeMap;
    for (auto&& elem : relevantCubeset.cubeMap()) {
        auto& cube = elem.first;
        REQUIRE(cube->tiledistance(0).genome() == 0);
        REQUIRE(cube->tiledistance(1).genome() == 1);
        REQUIRE(cube->tiledistance().size() == 2);
        auto hgSeq = idMap->querySequenceName(cube->tiledistance(0).sequence());
        auto mmSeq = idMap->querySequenceName(cube->tiledistance(1).sequence());
        if (orthology.at(hgSeq) == mmSeq) {
            orthologRelevantCubeMap[cube] = elem.second;
        } //else {
        //    std::cout << "[DEBUG] -- Discarding cube " << *cube << " with " << elem.second << std::endl;
        //}
    }

    //std::cout << "[DEBUG] -- idMap:" << *idMap << std::endl;
    for (auto&& elem : requiredRelevantCubeMap) {
        REQUIRE(orthologRelevantCubeMap.find(elem.first) != orthologRelevantCubeMap.end());
        REQUIRE(orthologRelevantCubeMap.at(elem.first) == elem.second);
    }
    REQUIRE(hopscotchEqual(requiredRelevantCubeMap, orthologRelevantCubeMap));
}
