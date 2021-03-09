#include <cstdlib>
#include <climits>

#include "catch2/catch.hpp"
#include "../SeedFinder.h"
#include "ConfigurationGenerator.h"
#include "LoadTestdata.h"

auto copyIDMap(std::shared_ptr<IdentifierMapping const> cIdMap) {
    auto idMap = std::make_shared<IdentifierMapping>();
    for (size_t gid = 0; gid < cIdMap->numGenomes(); ++gid) {
        idMap->queryGenomeID( cIdMap->queryGenomeName(gid) );
    }
    for (size_t sid = 0; sid < cIdMap->numSequences(); ++sid) {
        auto sequenceTuple = cIdMap->querySequenceTuple(sid);
        auto genome = cIdMap->queryGenomeName(sequenceTuple.gid);
        idMap->querySequenceID(sequenceTuple.sequence,
                               genome);
    }
    return idMap;
}

template<typename LinksetType>
auto reduceLinksetToResult(LinksetType const & linkset) {
    tsl::hopscotch_set<Link, LinkHash> reducedLinks{};
    for (auto&& elem : linkset.linkset()) {
        auto& link = elem.first;
        if ((link.first().genome() == 0) && link.second().genome() == 1) {
            Link l(std::vector<KmerOccurrence>{link.first(),
                                               link.second()},
                   link.span());
            reducedLinks.emplace(l);
        }
    }
    return reducedLinks;
}

auto loadResultJson(fs::path const & file, IdentifierMapping const & idMap) {
    tsl::hopscotch_set<Link, LinkHash> expectedLinks{};
    // load json
    std::ifstream is(file);
    json j;
    is >> j;
    // fill expected link
    for (auto&& elem : j) {
        auto occ0 = elem[0];
        auto occ1 = elem[1];
        Link l(std::vector<KmerOccurrence>{KmerOccurrence{(uint8_t)idMap.queryGenomeIDConst(occ0[2].get<std::string>()),
                                                          (uint32_t)idMap.querySequenceIDConst(occ0[3].get<std::string>(), occ0[2].get<std::string>()),
                                                          static_cast<size_t>(std::stoi(occ0[0].get<std::string>())),
                                                          static_cast<bool>(std::stoi(occ0[1].get<std::string>())), ""},
                                           KmerOccurrence{(uint8_t)idMap.queryGenomeIDConst(occ1[2].get<std::string>()),
                                                          (uint32_t)idMap.querySequenceIDConst(occ1[3].get<std::string>(), occ1[2].get<std::string>()),
                                                          static_cast<size_t>(std::stoi(occ1[0].get<std::string>())),
                                                          static_cast<bool>(std::stoi(occ1[1].get<std::string>())), ""}},
                                           static_cast<size_t>(std::stoi(occ0[4].get<std::string>())));
        expectedLinks.emplace(l);
    }
    return expectedLinks;
}

auto loadCubeJson(fs::path const & file, IdentifierMapping const & idMap, Configuration const & config) {
    tsl::hopscotch_map<std::shared_ptr<Cube const>,
                       tsl::hopscotch_set<LinkPtr, LinkPtrHash>,
                       CubePtrHash, CubePtrEqual> expectedCubeMap{};
    // load json
    std::ifstream is(file);
    json j;
    is >> j;
    // fill cubeMap
    for (auto&& elem : j.items()) {
        auto jcube = json::parse(elem.key());
        auto jlinks = elem.value();
        // create cube
        std::vector<Tiledistance> cubev;
        long long minDist = 0;
        for (auto&& occ : jcube) {
            cubev.emplace_back(idMap.queryGenomeIDConst(occ[0].get<std::string>()),
                               idMap.querySequenceIDConst(occ[1].get<std::string>(), occ[0].get<std::string>()),
                               static_cast<long long>(std::stoi(occ[2].get<std::string>())),
                               static_cast<bool>(std::stoi(occ[3].get<std::string>())));
            minDist = std::min(minDist, cubev.rbegin()->distance());
        }
        auto cubelink = Link();
        for (auto&& td : cubev) {
            cubelink.insertOccurrence(td.genome(), td.sequence(), (td.distance() + (-1*minDist)), td.reverse(), "A");
        }
        auto cube = std::make_shared<Cube const>(cubelink, config.tileSize());
        if (cube->tiledistance(0).distance() != 0) {
            std::cerr << "[ERROR] -- LoadTestdataMetagraph::expectedCubes -- cube '" << *cube << "' from '" << jcube << "'" << std::endl;
            throw std::runtime_error("[ERROR] -- LoadTestdataMetagraph::expectedCubes -- Cube creation corrupted");
        }
        std::sort(cubev.begin(), cubev.end());
        if (!(cubev == cube->tiledistance())) {
            std::cerr << "[ERROR] -- loadCubeJson -- cube vs. cubev: '" << cube << "' vs. '" << cubev << "'" << std::endl;
            throw std::runtime_error("[ERROR] -- loadCubeJson -- Cube creation corrupted");
        }
        // create links
        tsl::hopscotch_set<LinkPtr, LinkPtrHash> links;
        for (auto&& elem2 : jlinks.items()) {
            if (elem2.key()[0] != '[') { continue; }
            auto link = Link();
            auto jlink = json::parse(elem2.key());
            size_t span = 1;
            for (auto&& occ : jlink) {
                span = static_cast<size_t>(std::stoi(occ[4].get<std::string>()));
                link.insertOccurrence(idMap.queryGenomeIDConst(occ[2].get<std::string>()),
                                      idMap.querySequenceIDConst(occ[3].get<std::string>(), occ[2].get<std::string>()),
                                      static_cast<size_t>(std::stoi(occ[0].get<std::string>())),
                                      static_cast<bool>(std::stoi(occ[1].get<std::string>())),
                                      "A"); // kmer does not matter
            }
            link.extendSpanToRight(span - 1); // link initialized with span = 1
            links.emplace(LinkPtr(link));
        }
        // group links as links in expected cubes are grouped in python but not in SeedFinder::GH
        //std::vector<LinkPtr> sortedLinks(links.begin(), links.end());
        //links.clear();
        //std::sort(sortedLinks.begin(), sortedLinks.end());
        //groupOverlappingSeeds(sortedLinks);
        //links.insert(sortedLinks.begin(), sortedLinks.end());

        expectedCubeMap.emplace(cube, links);
    }
    return expectedCubeMap;
}



TEST_CASE("Seed Finder") {
    std::cout << "[INFO] -- [TEST CASE] -- Seed Finder" << std::endl;
    auto testdata = LoadTestdataMetagraph();
    auto hash = std::hash<std::string>{}(testdata.graphAnnotationFile());
    fs::path outfile{"/tmp/"+std::to_string(hash)+".json"};
    fs::path outfileCubes{"/tmp/"+std::to_string(hash)+".cubes.json"};
    auto cbuilder = ConfigBuilder();
    cbuilder.set(Allvsall{true},
                 Genome1{"genome0.fa"},
                 Genome2{"genome1.fa"},
                 Hasse(true),
                 MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{18},
                                                                              SpacedSeedMaskCollection::Span{18},
                                                                              SpacedSeedMaskCollection::SeedSetSize{1})},
                 MatchLimit{SIZE_MAX},
                 MetagraphInterfacePtr{testdata.metagraphInterface()},
                 OutputPath{outfile},
                 Redmask{true},
                 Thinning{5});

    SECTION("Contiguous, only seedFinding") {
        auto config = cbuilder.makeConfig();
        auto masks = config->maskCollection();
        REQUIRE(masks->size() == 1);
        REQUIRE(masks->weight() == masks->maxSpan());
        REQUIRE(masks->masksAsString() == std::vector<std::string>{"111111111111111111"});
        // Observe SeedFinder
        {
            auto observedSeedFinder = SeedFinder<TwoBitKmerDataShort>(config);
            observedSeedFinder.run();
        }
        // Explicit pipeline
        auto idMap = testdata.idMapGraph();
        auto observedResult = loadResultJson(outfile, *idMap);
        auto expectedSeedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, idMap);
        auto expectedLinkset = std::make_shared<Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>>(config, idMap);
        SeedKmerMap<TwoBitKmerDataShort> seedKmerMap{config, idMap};
        seedKmerMap.extractSeeds();
        expectedLinkset->createLinks(seedKmerMap);
        expectedLinkset->groupOverlappingLinks();
        auto expectedResult = reduceLinksetToResult(*expectedLinkset);
        REQUIRE(hopscotchSetEqual(expectedResult, observedResult));
    }

    SECTION("Four optimal spaced, with GH") {
        cbuilder.set(Allvsall{true},
                     CubeOutput{2},
                     CubeScoreThreshold{0},
                     Hasse(true),
                     MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{18},
                                                                                  SpacedSeedMaskCollection::SeedSetSize{4})},
                     MatchLimit{SIZE_MAX},
                     NThreads{1},
                     OptimalSeed{true},
                     PerformGeometricHashing{true},
                     TileSize{500});
        SECTION("Metagraph Input") {
            auto config = cbuilder.makeConfig();
            auto masks = config->maskCollection();
            REQUIRE(masks->size() == 4);
            REQUIRE(masks->weight() == 18);
            REQUIRE(masks->maxSpan() == 39);
            // Observe SeedFinder
            {
                auto observedSeedFinder = SeedFinder<TwoBitKmerDataShort>(config);
                observedSeedFinder.run();
            } // destruct to close outstreams
            auto idMap = testdata.idMapGraph();
            auto observedResult = loadResultJson(outfile, *idMap);
            auto observedCubes = loadCubeJson(outfileCubes, *idMap, *config);
            // SeedFinder pipeline
            auto expectedCubes = testdata.expectedCubesGraph(*idMap);
            REQUIRE(hopscotchEqual(expectedCubes, observedCubes));
            // Explicit pipeline rebuild
            auto expectedSeedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, idMap);
            auto expectedLinkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMap);
            SeedKmerMap<TwoBitKmerDataShort> seedKmerMap{config, idMap};
            seedKmerMap.extractSeeds();
            expectedLinkset->createLinks(seedKmerMap);
            Cubeset cubeset(expectedLinkset, testdata.expectedSequenceLengthsGraph());
            expectedLinkset->clear();
            for (auto&& elem : cubeset.scoreToCube()) {
                if (elem.first >= config->cubeScoreThreshold()) {
                    for (auto&& cube : elem.second) {
                        for (auto&& link : cubeset.cubeMap().at(cube)) {
                            expectedLinkset->createLinks(link.occurrence(), link.span());
                        }
                        if (cubeset.subcubeMap().find(cube) != cubeset.subcubeMap().end()) {
                            for (auto&& subcube : cubeset.subcubeMap().at(cube)) {
                                for (auto&& link : cubeset.cubeMap().at(subcube)) {
                                    expectedLinkset->createLinks(link.occurrence(), link.span());
                                }
                            }
                        }
                    }
                }
            }
            expectedLinkset->groupOverlappingLinks();
            auto expectedResult = reduceLinksetToResult(*expectedLinkset);
            REQUIRE(hopscotchSetEqual(expectedResult, observedResult));
        }
        SECTION("Fasta Input") {
            cbuilder.set(Genome1{testdata.configFasta()->genome1()},
                         Genome2{testdata.configFasta()->genome2()},
                         InputFiles{testdata.configFasta()->inputFiles()});
            auto config = cbuilder.makeConfig();
            auto masks = config->maskCollection();
            REQUIRE(masks->size() == 4);
            REQUIRE(masks->weight() == 18);
            REQUIRE(masks->maxSpan() == 39);
            // Observe SeedFinder
            {
                auto observedSeedFinder = SeedFinder<TwoBitKmerDataShort>(config);
                observedSeedFinder.run();
            } // destruct to close outstreams
            auto idMap = testdata.idMapFasta();
            auto observedResult = loadResultJson(outfile, *idMap);
            auto observedCubes = loadCubeJson(outfileCubes, *idMap, *config);
            // SeedFinder pipeline
            auto expectedCubes = testdata.expectedCubesFasta(*(testdata.idMapGraph()));
            REQUIRE(hopscotchEqual(expectedCubes, observedCubes));
            // Explicit pipeline rebuild
            auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, idMap);
            auto fastaCollection = std::make_shared<FastaCollection>(config);
            seedMap->extractSeeds(fastaCollection);
            auto expectedLinkset = std::make_shared<Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>>(config, idMap);
            expectedLinkset->createLinks(*seedMap);
            Cubeset cubeset(expectedLinkset, testdata.expectedSequenceLengthsExact());
            expectedLinkset->clear();
            for (auto&& elem : cubeset.scoreToCube()) {
                if (elem.first >= config->cubeScoreThreshold()) {
                    for (auto&& cube : elem.second) {
                        for (auto&& link : cubeset.cubeMap().at(cube)) {
                            expectedLinkset->createLinks(link.occurrence(), link.span());
                        }
                        if (cubeset.subcubeMap().find(cube) != cubeset.subcubeMap().end()) {
                            for (auto&& subcube : cubeset.subcubeMap().at(cube)) {
                                for (auto&& link : cubeset.cubeMap().at(subcube)) {
                                    expectedLinkset->createLinks(link.occurrence(), link.span());
                                }
                            }
                        }
                    }
                }
            }
            expectedLinkset->groupOverlappingLinks();
            auto expectedResult = reduceLinksetToResult(*expectedLinkset);
            REQUIRE(hopscotchSetEqual(expectedResult, observedResult));
        }
        SECTION("Fasta Input With Pre-GH") {
            cbuilder.set(Genome1{testdata.configFasta()->genome1()},
                         Genome2{testdata.configFasta()->genome2()},
                         InputFiles{testdata.configFasta()->inputFiles()},
                         PreHasse{true},
                         PreLinkThreshold{2},
                         PreMaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{20},
                                              SpacedSeedMaskCollection::SeedSetSize{2})},
                         PreOptimalSeed{true});
            auto config = cbuilder.makeConfig();
            auto masks = config->maskCollection();
            REQUIRE(masks->size() == 4);
            REQUIRE(masks->weight() == 18);
            REQUIRE(masks->maxSpan() == 39);
            REQUIRE(config->preMaskCollection()->size() == 2);
            REQUIRE(config->preMaskCollection()->weight() == 20);
            // Observe SeedFinder
            {
                auto observedSeedFinder = SeedFinder<TwoBitKmerDataShort>(config);
                observedSeedFinder.run();
            } // destruct to close outstreams
            auto idMap = testdata.idMapFasta();
            auto observedResult = loadResultJson(outfile, *idMap);
            auto observedCubes = loadCubeJson(outfileCubes, *idMap, *config);
            // TODO: actually test results, so far only test if it runs in principle
        }
        SECTION("Fasta Input With Pre-Filter, M3") {
            cbuilder.set(Genome1{testdata.configFasta()->genome1()},
                         Genome2{testdata.configFasta()->genome2()},
                         InputFiles{testdata.configFasta()->inputFiles()},
                         PerformGeometricHashing{false},
                         PreLinkThreshold{2},
                         PreMaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{20},
                                              SpacedSeedMaskCollection::SeedSetSize{2})},
                         PreOptimalSeed{true});
            auto config = cbuilder.makeConfig();
            auto masks = config->maskCollection();
            REQUIRE(masks->size() == 4);
            REQUIRE(masks->weight() == 18);
            REQUIRE(masks->maxSpan() == 39);
            REQUIRE(config->preMaskCollection()->size() == 2);
            REQUIRE(config->preMaskCollection()->weight() == 20);
            // Observe SeedFinder
            {
                auto observedSeedFinder = SeedFinder<TwoBitKmerDataShort>(config);
                observedSeedFinder.run();
            } // destruct to close outstreams
            auto idMap = testdata.idMapFasta();
            auto observedResult = loadResultJson(outfile, *idMap);
            auto observedCubes = loadCubeJson(outfileCubes, *idMap, *config);
            // TODO: actually test results, so far only test if it runs in principle
        }
    }

    SECTION("Four optimal spaced, with GH, 1-vs-all") {
        auto hash = std::hash<std::string>{}(testdata.graphAnnotationFile());
        fs::path outfile{"/tmp/"+std::to_string(hash)+".json"};
        fs::path outfileCubes{"/tmp/"+std::to_string(hash)+".cubes.json"};
        cbuilder.set(Allvsall{false},
                     CubeOutput{2},
                     CubeScoreThreshold{0},
                     MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{18},
                                                                                  SpacedSeedMaskCollection::SeedSetSize{4})},
                     MatchLimit{SIZE_MAX},
                     MetagraphInterfacePtr{testdata.metagraphInterface()},
                     OptimalSeed{true},
                     OutputPath{outfile},
                     PerformGeometricHashing{true},
                     TileSize{500});
        auto config = cbuilder.makeConfig();
        auto masks = config->maskCollection();
        REQUIRE(masks->size() == 4);
        REQUIRE(masks->weight() == 18);
        REQUIRE(masks->maxSpan() == 39);
        // Observe SeedFinder
        {
            auto observedSeedFinder = SeedFinder<TwoBitKmerDataShort>(config);
            observedSeedFinder.run();
        } // destruct to close outstreams
        auto idMap = *(testdata.idMapGraph());
        auto observedCubes = loadCubeJson(outfileCubes, idMap, *config);

        // Explicit pipeline
        auto expectedCubes = testdata.expectedCubesGraph(idMap);
        REQUIRE(expectedCubes == observedCubes);
    }

    SECTION("Four optimal spaced, with GH, batch mode") {
        auto hash = std::hash<std::string>{}(testdata.graphAnnotationFile());
        fs::path outfile{"/tmp/"+std::to_string(hash)+".json"};
        fs::path outfileCubes{"/tmp/"+std::to_string(hash)+".cubes.json"};
        cbuilder.set(Allvsall{true},
                     Batchsize{2},
                     CubeOutput{2},
                     CubeScoreThreshold{0},
                     Genome1{testdata.configFasta()->genome1()},
                     Genome2{testdata.configFasta()->genome2()},
                     Hasse(true),
                     InputFiles{testdata.configFasta()->inputFiles()},
                     MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{18},
                                                                                  SpacedSeedMaskCollection::SeedSetSize{4})},
                     MatchLimit{SIZE_MAX},
                     OptimalSeed{true},
                     OutputPath{outfile},
                     PerformGeometricHashing{true},
                     TileSize{500});
        auto config = cbuilder.makeConfig();
        auto masks = config->maskCollection();
        REQUIRE(masks->size() == 4);
        REQUIRE(masks->weight() == 18);
        REQUIRE(masks->maxSpan() == 39);
        // Observe SeedFinder
        auto idMap = testdata.idMapFasta();
        std::vector<std::vector<size_t>> observedBatches;
        {
            auto observedSeedFinder = SeedFinder<TwoBitKmerDataShort>(config);
            observedSeedFinder.run();
            observedBatches = observedSeedFinder.fastaBatches(*idMap);
        } // destruct to close outstreams
        auto observedCubes = loadCubeJson(outfileCubes, *idMap, *config);
        // check batches
        auto allPairings = [](std::vector<size_t> sids, IdentifierMapping const & idMap) {
            std::set<std::vector<size_t>> pairings;
            std::vector<std::vector<size_t>> seqIDs(idMap.numGenomes());
            for (auto&& sid : sids) {
                seqIDs.at(idMap.queryGenomeIDFromSequence(sid)).emplace_back(sid);
            }
            for (auto&& vec : seqIDs) {
                if (pairings.size() == 0) {
                    for (auto&& sid : vec) { pairings.insert(std::vector<size_t>{sid}); }
                } else {
                    std::set<std::vector<size_t>> newPairings;
                    for (auto&& pair : pairings) {
                        for (auto&& sid : vec) {
                            auto newpair = pair;
                            newpair.emplace_back(sid);
                            newPairings.insert(newpair);
                        }
                    }
                    pairings = newPairings;
                }
            }
            return pairings;
        };
        std::vector<size_t> allSids;
        for (auto&& sid : idMap->sequenceIDToTuple()) { allSids.emplace_back(sid.first); }
        auto expectedPairings = allPairings(allSids, *idMap);
        std::set<std::vector<size_t>> observedPairings;
        for (auto&& batch : observedBatches) {
            auto batchPairings = allPairings(batch, *idMap);
            for (auto&& pairing : batchPairings) {
                REQUIRE(expectedPairings.find(pairing) != expectedPairings.end());
                REQUIRE(observedPairings.find(pairing) == observedPairings.end());
                observedPairings.insert(pairing);
            }
        }
        REQUIRE(expectedPairings == observedPairings);
        // Explicit pipeline
        auto expectedCubes = testdata.expectedCubesFasta(*(testdata.idMapGraph()));
        //REQUIRE(expectedCubes == observedCubes);
        for (auto&& linkset : observedCubes) {
            REQUIRE(hopscotchSetEqual(linkset.second, expectedCubes.at(linkset.first)));
        }
        for (auto&& linkset : expectedCubes) {
            REQUIRE(hopscotchSetEqual(linkset.second, observedCubes.at(linkset.first)));
        }
        REQUIRE(hopscotchEqual(observedCubes, expectedCubes));
    }
}
