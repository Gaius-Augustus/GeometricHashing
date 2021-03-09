#include <experimental/filesystem>
#include <thread>
#include <fstream>
#include <memory>

#include "catch2/catch.hpp"
#include "../IdentifierMapping.h"
#include "../Linkset.h"
#include "../SeedMap.h"
#include "fillLinksetTestdata.h"
#include "GeometricHashingTestdata.h"
#include "LoadTestdata.h"



TEST_CASE("Test Cartesian Product Helper Function") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Cartesian Product Helper Function" << std::endl;
    std::vector<std::vector<size_t>> input{std::vector<size_t>{1,2,3},
                                           std::vector<size_t>{10,20,30},
                                           std::vector<size_t>{100,200,300}};
    for (size_t i = 0; i < 3*3*3; ++i) {
        if (i == 0) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,10,100}); }
        if (i == 1) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,10,100}); }
        if (i == 2) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,10,100}); }
        if (i == 3) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,20,100}); }
        if (i == 4) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,20,100}); }
        if (i == 5) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,20,100}); }
        if (i == 6) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,30,100}); }
        if (i == 7) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,30,100}); }
        if (i == 8) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,30,100}); }
        if (i == 9) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,10,200}); }
        if (i == 10) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,10,200}); }
        if (i == 11) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,10,200}); }
        if (i == 12) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,20,200}); }
        if (i == 13) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,20,200}); }
        if (i == 14) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,20,200}); }
        if (i == 15) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,30,200}); }
        if (i == 16) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,30,200}); }
        if (i == 17) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,30,200}); }
        if (i == 18) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,10,300}); }
        if (i == 19) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,10,300}); }
        if (i == 20) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,10,300}); }
        if (i == 21) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,20,300}); }
        if (i == 22) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,20,300}); }
        if (i == 23) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,20,300}); }
        if (i == 24) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,30,300}); }
        if (i == 25) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{2,30,300}); }
        if (i == 26) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{3,30,300}); }
        if (i == 27) { REQUIRE(cartesianProductByID(i, input) == std::vector<size_t>{1,10,100}); } // overflow
    };
}



TEST_CASE("Test Linkset on normal Testdata") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Linkset on normal Testdata" << std::endl;
    LoadTestdata data;
    auto fastaCollection = data.fastaCollection();
    auto cbuilder = data.cbuilder();
    cbuilder.set(Hasse(true),
                 MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{20},
                                                                              SpacedSeedMaskCollection::Span{20},
                                                                              SpacedSeedMaskCollection::SeedSetSize{1})},
                 MatchLimit{SIZE_MAX},
                 MatchLimitDiscardSeeds{false},
                 OccurrencePerGenomeMax{1},
                 OccurrencePerGenomeMin{0},
                 OccurrencePerSequenceMax{SIZE_MAX},
                 Thinning{1});
    auto config = cbuilder.makeConfig();
    auto kmerMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, data.idMap());
    kmerMap->extractSeeds(fastaCollection);

    SECTION("Test merging") {
        auto idMap = std::make_shared<IdentifierMapping>(*(kmerMap->idMap()));
        // make sure enough ids are known
        idMap->queryGenomeID("gen0");
        idMap->queryGenomeID("gen1");
        idMap->queryGenomeID("gen2");
        idMap->querySequenceID("s0", "gen0");
        idMap->querySequenceID("s1", "gen0");
        idMap->querySequenceID("s2", "gen0");
        idMap->querySequenceID("s0", "gen1");
        idMap->querySequenceID("s1", "gen1");
        idMap->querySequenceID("s2", "gen1");
        idMap->querySequenceID("s0", "gen2");
        idMap->querySequenceID("s1", "gen2");
        idMap->querySequenceID("s2", "gen2");
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> l1{config, idMap};
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> l2{config, idMap};
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> l3{config, idMap};
        for (auto&& occurrences : {std::vector<KmerOccurrence>{KmerOccurrence{0, 0, 100, false, "A"},
                                                               KmerOccurrence{1, 3, 100, false, "A"},
                                                               KmerOccurrence{2, 6, 100, false, "A"}},
                                   std::vector<KmerOccurrence>{KmerOccurrence{0, 1, 100, false, "A"},
                                                               KmerOccurrence{1, 4, 100, false, "A"},
                                                               KmerOccurrence{2, 7, 100, false, "A"}},
                                   std::vector<KmerOccurrence>{KmerOccurrence{0, 2, 100, false, "A"},
                                                               KmerOccurrence{1, 5, 100, false, "A"},
                                                               KmerOccurrence{2, 8, 100, false, "A"}}}) {
            l1.createLinks(occurrences, config->span());
            l3.createLinks(occurrences, config->span());
        }
        for (auto&& occurrences : {std::vector<KmerOccurrence>{KmerOccurrence{0, 0, 200, false, "A"},
                                                               KmerOccurrence{1, 3, 500, false, "A"},
                                                               KmerOccurrence{2, 6, 100, false, "A"}},
                                   std::vector<KmerOccurrence>{KmerOccurrence{0, 1, 1000, false, "A"},
                                                               KmerOccurrence{1, 4, 700, false, "A"},
                                                               KmerOccurrence{2, 8, 200, false, "A"}}}) {
            l2.createLinks(occurrences, config->span());
            l3.createLinks(occurrences, config->span());
        }
        l1.merge(l2);
        REQUIRE(l1 == l3);
    }

    SECTION("Test with restriction max_occurrences_per_genome = 1") {
//std::cout << std::endl << std::endl << std::endl << std::endl << "[DEBUG] 000" << std::endl << std::endl << std::endl << std::endl << std::endl;
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset(config, kmerMap->idMap());
        auto idMap = *(linkset.idMapping()); //constIDMapRef; // true copy so that original mapping remains unchanged
        linkset.createLinks(*kmerMap);
        auto& lSet = linkset.linkset();

        //REQUIRE(linkset.numDiscardedKmers() == 113);
        REQUIRE(linkset.size() == 166);
        REQUIRE(linkset.size() == lSet.size());

        // create links that are expected and check if (only) they occur
        typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType expectedLinkset;
        fillLinksetTestdata(expectedLinkset, idMap);

        // check if only expected links occur
        for (auto&& expected : expectedLinkset) {
            REQUIRE(lSet.find(expected.first) != lSet.end());
        }
        // check if all expected links occur
        for (auto&& observed : lSet) {
            auto obsLink = observed.first;
            REQUIRE(expectedLinkset.find(obsLink) != expectedLinkset.end());
        }
    }
}



TEST_CASE("Test Linkset on geometricHashing Testdata") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Linkset on geometricHashing Testdata" << std::endl;
    auto data = LoadTestdata(LoadTestdata::geometricHashingData);
    auto fastaCollection = data.fastaCollection();
    auto cbuilder = data.cbuilder();
    cbuilder.set(Hasse(true),
                 MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection>(SpacedSeedMaskCollection::Weight{5},
                                                                              SpacedSeedMaskCollection::Span{5},
                                                                              SpacedSeedMaskCollection::SeedSetSize{1})},
                 MatchLimit{SIZE_MAX},
                 MatchLimitDiscardSeeds{false},
                 OccurrencePerGenomeMax{1},
                 OccurrencePerGenomeMin{1},
                 Thinning{1});
    auto config = cbuilder.makeConfig();
    auto kmerMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, data.idMap());
    kmerMap->extractSeeds(fastaCollection);

    SECTION("Test with restriction max_occurrences_per_genome = 1") {
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset(config, kmerMap->idMap());
        linkset.createLinks(*kmerMap);

        auto idMap = *(linkset.idMapping()); //constIDMapRef; // true copy so that original mapping remains unchanged
        auto& lSet = linkset.linkset();
        REQUIRE(linkset.numDiscardedKmers() == 1);
        REQUIRE(linkset.size() == 15);  // as now exact positions are used in Links, TTGGC and TATGC create different Links
        REQUIRE(linkset.size() == lSet.size());
        //REQUIRE(idMap.refGenome() == kTestSpecies[0]);
        REQUIRE(idMap.queryGenomeID(kTestSpecies[0]) == static_cast<uint16_t>(0));
        REQUIRE(idMap.queryGenomeName(0) == kTestSpecies[0]);
        REQUIRE(idMap.querySequenceIDConst(kTestSequences[0], kTestSpecies[0]) == static_cast<uint16_t>(0));
        REQUIRE(idMap.querySequenceName(0) == (kTestSequences[0]));

        // create links that are expected and check if (only) they occur
        ManualSets manSet(idMap);
        auto& expectedLinkset = manSet.expectedLinkset();
        // check if only expected links occur
        for (auto&& expected : expectedLinkset) {
            REQUIRE(lSet.find(expected) != lSet.end());
        }
        // check if all expected links occur
        for (auto&& observed : lSet) {
            auto obsLink = observed.first;
            REQUIRE(expectedLinkset.find(obsLink) != expectedLinkset.end());
        }

        SECTION("Test Link observation counts") {
            auto duplicateLinkNc = Link();
            // Link 00-11 (pos. 0, 0)
            duplicateLinkNc.insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]),
                                             idMap.querySequenceID(kTestSequences[0], kTestSpecies[0]),
                                             0, false, "ACGT");
            duplicateLinkNc.insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]),
                                             idMap.querySequenceID(kTestSequences[1], kTestSpecies[1]),
                                             0, false, "ACGT");
            duplicateLinkNc.extendSpanToRight(4);
            Link const duplicateLink = duplicateLinkNc;    // comparisons overloaded with <Link const>
            REQUIRE(linkset.linkCount(duplicateLink) == 1); // from data
            linkset.addLink(duplicateLink);
            REQUIRE(linkset.linkCount(duplicateLink) == 2); // 'manual' increase
            // check other Link counts
            for (auto&& observed : lSet) {
                auto obsLink = observed.first;
                auto timesObserved = observed.second;
                if (obsLink == duplicateLink) { continue; }
                REQUIRE(timesObserved == 1);
                REQUIRE(linkset.linkCount(obsLink) == 1);
            }
        }
    }

    SECTION("Test with restriction min_occurrences_per_genome = 2") {
        cbuilder.set(OccurrencePerGenomeMax{SIZE_MAX},
                     OccurrencePerGenomeMin{2});
        auto config2 = cbuilder.makeConfig();
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset(config2, kmerMap->idMap());
        linkset.createLinks(*kmerMap);

        auto idMap = *(linkset.idMapping());//constIDMapRef; // true copy so that original mapping remains unchanged
        auto& lSet = linkset.linkset();

        REQUIRE(linkset.numDiscardedKmers() == 15); // 00 - 10 from two different k-mers
        REQUIRE(linkset.size() == 4);
        for (auto&& linkToCount : lSet) {
            auto count = linkToCount.second;
            REQUIRE(count == static_cast<size_t>(1));
        }

        // create links that are expected and check if (only) they occur
        auto l1 = Link();
        auto l2 = Link();
        auto l3 = Link();
        auto l4 = Link();
        // exact positions stored
        l1.insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]), idMap.querySequenceID(kTestSequences[0], kTestSpecies[0]), 90, false, "ACGT");
        l1.insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]), idMap.querySequenceID(kTestSequences[1], kTestSpecies[1]), 42, false, "ACGT");
        l2.insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]), idMap.querySequenceID(kTestSequences[0], kTestSpecies[0]), 105, false, "ACGT");
        l2.insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]), idMap.querySequenceID(kTestSequences[1], kTestSpecies[1]), 105, false, "ACGT");
        l3.insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]), idMap.querySequenceID(kTestSequences[0], kTestSpecies[0]), 90, false, "ACGT");
        l3.insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]), idMap.querySequenceID(kTestSequences[1], kTestSpecies[1]), 105, false, "ACGT");
        l4.insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]), idMap.querySequenceID(kTestSequences[0], kTestSpecies[0]), 105, false, "ACGT");
        l4.insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]), idMap.querySequenceID(kTestSequences[1], kTestSpecies[1]), 42, false, "ACGT");
        l1.extendSpanToRight(4);
        l2.extendSpanToRight(4);
        l3.extendSpanToRight(4);
        l4.extendSpanToRight(4);
        std::unordered_set<Link, LinkHash> expectedLinkset{l1, l2, l3, l4};
        // check if only expected links occur
        for (auto&& expected : expectedLinkset) {
            REQUIRE(lSet.find(expected) != lSet.end());
        }
        // check if all expected links occur
        for (auto&& observed : lSet) {
            auto obsLink = observed.first;
            REQUIRE(expectedLinkset.find(obsLink) != expectedLinkset.end());
        }
    }
}



TEST_CASE("Test Expected Links with Metagraph") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Expected Links with Metagraph" << std::endl;
    auto data = LoadTestdataMetagraph();
    auto linksetFasta = data.runLinksetCreationFasta();
    auto linksetGraph = data.runLinksetCreationGraph();
    auto expectedLinksFasta = data.expectedLinksFasta(*data.idMapGraph());
    auto expectedLinksGraph = data.expectedLinksGraph(*data.idMapGraph());
    REQUIRE(hopscotchEqual(expectedLinksFasta, linksetFasta->linkset()));
std::cout << *(data.idMapGraph()) << std::endl;
    REQUIRE(hopscotchEqual(expectedLinksGraph, linksetGraph->linkset()));
}




TEST_CASE("Match Limit Filter") {
    std::cout << "[INFO] -- [TEST CASE] -- Match Limit Filter" << std::endl;
    auto idMap = std::make_shared<IdentifierMapping>("genome0");
    idMap->queryGenomeID("genome1");
    idMap->querySequenceID("seq0", "genome0");
    idMap->querySequenceID("seq1", "genome1");
    SECTION("Two Genomes") {
        auto cbuilder = ConfigBuilder();
        cbuilder.set(CreateAllMatches{false},
                     Genome1{"genome0"},
                     Genome2{"genome1"},
                     Hasse(true),
                     MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{4},
                                                                                        SpacedSeedMaskCollection::Span{4},
                                                                                        SpacedSeedMaskCollection::SeedSetSize{1})},
                     MatchLimit{2},
                     MatchLimitDiscardSeeds{false});
        auto config = cbuilder.makeConfig();
        cbuilder.set(CreateAllMatches{true},
                     MatchLimitDiscardSeeds{false});
        auto config2 = cbuilder.makeConfig();
        cbuilder.set(CreateAllMatches{false},
                     MatchLimitDiscardSeeds{true});
        auto config3 = cbuilder.makeConfig();
        // do not create all
        auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, idMap);
        // create all, should not make a difference in this setting
        auto seedMap2 = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config2, idMap);
        // discard filtered
        auto seedMap3 = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config3, idMap);
        auto o11 = KmerOccurrence(0, 0, 10, false, "ACGT");
        auto o12 = KmerOccurrence(0, 0, 20, false, "ACGT");
        auto o13 = KmerOccurrence(0, 0, 30, false, "ACGT");
        auto o14 = KmerOccurrence(0, 0, 40, false, "ACGT");
        auto o21 = KmerOccurrence(1, 1, 15, false, "ACGT");
        auto o22 = KmerOccurrence(1, 1, 25, false, "ACGT");
        auto o23 = KmerOccurrence(1, 1, 35, false, "ACGT");
        auto o24 = KmerOccurrence(1, 1, 45, false, "ACGT");
        seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // six matches, filtered to two
        seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
        seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o13, 0);
        seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
        seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
        seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o14, 0); // two matches
        seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
        seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o24, 0);
        // should be equal
        seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // six matches, filtered to two
        seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
        seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o13, 0);
        seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
        seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
        seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o14, 0); // two matches
        seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
        seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o24, 0);
        // should only contain last two matches
        seedMap3->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // six matches, discarded
        seedMap3->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
        seedMap3->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o13, 0);
        seedMap3->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
        seedMap3->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
        seedMap3->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o14, 0); // two matches
        seedMap3->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
        seedMap3->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o24, 0);

        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset{config, idMap};
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset2{config2, idMap};
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset3{config3, idMap};

        linkset.createLinks(*seedMap);
        linkset2.createLinks(*seedMap2);
        linkset3.createLinks(*seedMap3);

        REQUIRE(linkset.size() == 4);
        REQUIRE(linkset2.size() == 4);
        REQUIRE(linkset3.size() == 2);

        auto allPossible = std::set<Link>{Link(std::vector<KmerOccurrence>{o11, o21}, config->span()),
                                          Link(std::vector<KmerOccurrence>{o11, o22}, config->span()),
                                          Link(std::vector<KmerOccurrence>{o12, o21}, config->span()),
                                          Link(std::vector<KmerOccurrence>{o12, o22}, config->span()),
                                          Link(std::vector<KmerOccurrence>{o13, o21}, config->span()),
                                          Link(std::vector<KmerOccurrence>{o13, o22}, config->span()),
                                          Link(std::vector<KmerOccurrence>{o14, o23}, config->span()),
                                          Link(std::vector<KmerOccurrence>{o14, o24}, config->span())};
        REQUIRE(linkset.linkset().find(Link(std::vector<KmerOccurrence>{o14, o23}, config->span())) != linkset.linkset().end()); // find two unfiltered matches
        REQUIRE(linkset.linkset().find(Link(std::vector<KmerOccurrence>{o14, o24}, config->span())) != linkset.linkset().end());
        for (auto&& elem : linkset.linkset()) {
            auto& match = elem.first;
            REQUIRE(allPossible.find(match) != allPossible.end());  // remaining matches must be from the set of possible matches
        }
        REQUIRE(linkset2.linkset().find(Link(std::vector<KmerOccurrence>{o14, o23}, config->span())) != linkset2.linkset().end());
        REQUIRE(linkset2.linkset().find(Link(std::vector<KmerOccurrence>{o14, o24}, config->span())) != linkset2.linkset().end());
        for (auto&& elem : linkset2.linkset()) {
            auto& match = elem.first;
            REQUIRE(allPossible.find(match) != allPossible.end());
        }
        REQUIRE(linkset3.linkset().find(Link(std::vector<KmerOccurrence>{o14, o23}, config->span())) != linkset3.linkset().end());
        REQUIRE(linkset3.linkset().find(Link(std::vector<KmerOccurrence>{o14, o24}, config->span())) != linkset3.linkset().end());
    }

    SECTION("More Genomes") {
        idMap->queryGenomeID("genome2");
        idMap->queryGenomeID("genome3");
        idMap->querySequenceID("seq2", "genome2");
        idMap->querySequenceID("seq3", "genome3");

        auto o11 = KmerOccurrence(0, 0, 10, false, "ACGT");
        auto o12 = KmerOccurrence(0, 0, 20, false, "ACGT");
        auto o21 = KmerOccurrence(1, 1, 15, false, "ACGT");
        auto o22 = KmerOccurrence(1, 1, 25, false, "ACGT");
        auto o31 = KmerOccurrence(2, 2, 11, false, "ACGT");
        auto o32 = KmerOccurrence(2, 2, 21, false, "ACGT");
        auto o41 = KmerOccurrence(3, 3, 16, false, "ACGT");
        auto o42 = KmerOccurrence(3, 3, 26, false, "ACGT");

        auto o13 = KmerOccurrence(0, 0, 30, false, "ACGA");
        auto o23 = KmerOccurrence(1, 1, 35, false, "ACGA");
        auto o43 = KmerOccurrence(3, 3, 36, false, "ACGA");

        auto cbuilder = ConfigBuilder();
        cbuilder.set(CreateAllMatches{true}, // does not matter anymore
                     Genome1{"genome0"},
                     Genome2{"genome1"},
                     Hasse(true),
                     MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{4},
                                                                                        SpacedSeedMaskCollection::Span{4},
                                                                                        SpacedSeedMaskCollection::SeedSetSize{1})},
                     MatchLimit{3},
                     MatchLimitDiscardSeeds{false});

        SECTION("Create All") {
            auto config = cbuilder.makeConfig();
            // sample
            auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, idMap);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // sample 3 from here
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o31, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o32, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o41, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o42, 0);

            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o13, 0); // make all (one link)
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o43, 0);
            Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset{config, idMap};
            linkset.createLinks(*seedMap);

            // discardauto cbuilder = ConfigBuilder();
            cbuilder.set(CreateAllMatches{true},
                         MatchLimitDiscardSeeds{true});
            auto config2 = cbuilder.makeConfig();
            auto seedMap2 = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config2, idMap);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // discard
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o31, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o32, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o41, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o42, 0);

            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o13, 0); // make all (one link)
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o43, 0);
            Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset2{config2, idMap};
            linkset2.createLinks(*seedMap2);

            auto allPossible = std::set<Link>{Link(std::vector<KmerOccurrence>{o11, o21, o31, o41}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o11, o21, o31, o42}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o11, o21, o32, o41}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o11, o21, o32, o42}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o11, o22, o31, o41}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o11, o22, o31, o42}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o11, o22, o32, o41}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o11, o22, o32, o42}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o12, o21, o31, o41}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o12, o21, o31, o42}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o12, o21, o32, o41}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o12, o21, o32, o42}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o12, o22, o31, o41}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o12, o22, o31, o42}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o12, o22, o32, o41}, config->span()),
                                              Link(std::vector<KmerOccurrence>{o12, o22, o32, o42}, config->span()),

                                              Link(std::vector<KmerOccurrence>{o13, o23, o43}, config->span())};

            REQUIRE(linkset.size() == 4);
            REQUIRE(linkset.linkset().find(Link(std::vector<KmerOccurrence>{o13, o23, o43}, config->span())) != linkset.linkset().end()); // find unfiltered matches
            for (auto&& elem : linkset.linkset()) {
                auto& match = elem.first;
                REQUIRE(allPossible.find(match) != allPossible.end());  // remaining matches must be from the set of possible matches
            }
            REQUIRE(linkset2.size() == 1);
            REQUIRE(linkset2.linkset().find(Link(std::vector<KmerOccurrence>{o13, o23, o43}, config->span())) != linkset2.linkset().end()); // find unfiltered matches
        }
    }
}



TEST_CASE("Test Grouping of Overlapping Seeds") {
    std::cout << "[INFO] -- [TEST CASE] -- Grouping of Overlapping Seeds" << std::endl;
    SECTION("Pairwise") {
        tsl::hopscotch_set<Link, LinkHash> rawSet;
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,1,10,false,"AAAAA"}}, 5);
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,4,false,"AAAAA"}, KmerOccurrence{1,1,14,false,"AAAAA"}}, 5); // group with above
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,4,false,"AAAAA"}, KmerOccurrence{1,1,15,false,"AAAAA"}}, 5); // different diagonal, do not group
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,8,false,"AAAAA"}, KmerOccurrence{1,1,18,false,"AAAAA"}}, 5); // group with first two
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,13,false,"AAAAA"}, KmerOccurrence{1,1,23,false,"AAAAA"}}, 5); // not overlapping, do not group with first group
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,18,false,"AAAAA"}, KmerOccurrence{1,1,28,false,"AAAAA"}}, 5); // not overlapping, do not group with any above
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,20,false,"AAAAA"}, KmerOccurrence{1,1,30,false,"AAAAA"}}, 5); // overlapping, group

        tsl::hopscotch_set<Link, LinkHash> expectedSet;
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,1,10,false,"AAAAA"}}, 13); // last added at 8, span 5 -> 12 is last pos, thus span = 13
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,4,false,"AAAAA"}, KmerOccurrence{1,1,15,false,"AAAAA"}}, 5); // different diagonal, do not group
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,13,false,"AAAAA"}, KmerOccurrence{1,1,23,false,"AAAAA"}}, 5); // not overlapping, do not group with first group
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,18,false,"AAAAA"}, KmerOccurrence{1,1,28,false,"AAAAA"}}, 7); // second group, added up to pos 24, thus span = 7

        auto v = std::vector<Link>(rawSet.begin(), rawSet.end());
        std::sort(v.begin(), v.end());
        groupOverlappingSeeds(v);
        tsl::hopscotch_set<Link, LinkHash> observedSet{v.begin(), v.end()};
        REQUIRE(expectedSet == observedSet);
    }
    SECTION("Multiple") {
        tsl::hopscotch_set<Link, LinkHash> rawSet;
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,1,10,false,"AAAAA"}, KmerOccurrence{2,2,20,false,"AAAAA"}, KmerOccurrence{3,3,30,false,"AAAAA"}}, 5);
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,4,false,"AAAAA"}, KmerOccurrence{1,1,14,false,"AAAAA"}, KmerOccurrence{2,2,24,false,"AAAAA"}, KmerOccurrence{3,3,34,false,"AAAAA"}}, 5);  // group with above
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,8,false,"AAAAA"}, KmerOccurrence{1,1,18,false,"AAAAA"}, KmerOccurrence{2,2,28,false,"AAAAA"}, KmerOccurrence{3,3,38,false,"AAAAA"}}, 5);  // group with above
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,2,10,false,"AAAAA"}, KmerOccurrence{2,2,20,false,"AAAAA"}, KmerOccurrence{3,3,30,false,"AAAAA"}}, 5);  // different seq
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,1,10,false,"AAAAA"}, KmerOccurrence{2,2,20,false,"AAAAA"}, KmerOccurrence{4,4,30,false,"AAAAA"}}, 5);  // different genome
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,1,10,false,"AAAAA"}, KmerOccurrence{2,2,20,true,"AAAAA"}, KmerOccurrence{3,3,30,false,"AAAAA"}}, 5);   // different strand
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,12,false,"AAAAA"}, KmerOccurrence{1,1,23,false,"AAAAA"}, KmerOccurrence{2,2,32,false,"AAAAA"}, KmerOccurrence{3,3,42,false,"AAAAA"}}, 5); // different diagonal
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,12,false,"AAAAA"}, KmerOccurrence{1,1,22,false,"AAAAA"}, KmerOccurrence{2,2,32,false,"AAAAA"}, KmerOccurrence{3,3,42,false,"AAAAA"}}, 5); // group with first group
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,17,false,"AAAAA"}, KmerOccurrence{1,1,27,false,"AAAAA"}, KmerOccurrence{2,2,37,false,"AAAAA"}, KmerOccurrence{3,3,47,false,"AAAAA"}}, 5); // not overlapping
        rawSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,19,false,"AAAAA"}, KmerOccurrence{1,1,29,false,"AAAAA"}, KmerOccurrence{2,2,39,false,"AAAAA"}, KmerOccurrence{3,3,49,false,"AAAAA"}}, 5); // group with above

        tsl::hopscotch_set<Link, LinkHash> expectedSet;
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,1,10,false,"AAAAA"}, KmerOccurrence{2,2,20,false,"AAAAA"}, KmerOccurrence{3,3,30,false,"AAAAA"}}, 17);
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,2,10,false,"AAAAA"}, KmerOccurrence{2,2,20,false,"AAAAA"}, KmerOccurrence{3,3,30,false,"AAAAA"}}, 5);  // different seq
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,1,10,false,"AAAAA"}, KmerOccurrence{2,2,20,false,"AAAAA"}, KmerOccurrence{4,4,30,false,"AAAAA"}}, 5);  // different genome
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAAA"}, KmerOccurrence{1,1,10,false,"AAAAA"}, KmerOccurrence{2,2,20,true,"AAAAA"}, KmerOccurrence{3,3,30,false,"AAAAA"}}, 5);   // different strand
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,12,false,"AAAAA"}, KmerOccurrence{1,1,23,false,"AAAAA"}, KmerOccurrence{2,2,32,false,"AAAAA"}, KmerOccurrence{3,3,42,false,"AAAAA"}}, 5); // different diagonal
        expectedSet.emplace(std::vector<KmerOccurrence>{KmerOccurrence{0,0,17,false,"AAAAA"}, KmerOccurrence{1,1,27,false,"AAAAA"}, KmerOccurrence{2,2,37,false,"AAAAA"}, KmerOccurrence{3,3,47,false,"AAAAA"}}, 7); // not overlapping

        auto v = std::vector<Link>(rawSet.begin(), rawSet.end());
        std::sort(v.begin(), v.end());
        groupOverlappingSeeds(v);
        tsl::hopscotch_set<Link, LinkHash> observedSet{v.begin(), v.end()};
        REQUIRE(expectedSet == observedSet);
    }
}



TEST_CASE("Test Relevant Link Creation") {
    std::cout << "[INFO] -- [TEST CASE] -- Test Relevant Link Creation" << std::endl;
    auto idMap = std::make_shared<IdentifierMapping>("genome0");
    idMap->queryGenomeID("genome1");
    idMap->queryGenomeID("genome2");
    idMap->queryGenomeID("genome3");
    idMap->querySequenceID("seq0", "genome0");
    idMap->querySequenceID("seq1", "genome1");
    idMap->querySequenceID("seq2", "genome2");
    idMap->querySequenceID("seq3", "genome3");

    ConfigBuilder configBuilder{};
    configBuilder.set(Genome1{"genome0"},
                      Genome2{"genome1"},
                      Hasse(true),
                      NThreads{3},
                      MatchLimit{ULLONG_MAX},
                      OccurrencePerSequenceMax{ULLONG_MAX},
                      OptimalSeed{true},
                      PerformGeometricHashing{true},
                      PreHasse{true},
                      PreLinkThreshold{0},
                      PreOptimalSeed{true},
                      TileSize{50});

    tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> relevantCubes;
    relevantCubes.insert(std::make_shared<Cube>(std::vector<Tiledistance>{Tiledistance{0,0,0,false}, Tiledistance{1,1,0,false}}));
    relevantCubes.insert(std::make_shared<Cube>(std::vector<Tiledistance>{Tiledistance{0,0,0,false}, Tiledistance{1,1,-4,false}}));
    relevantCubes.insert(std::make_shared<Cube>(std::vector<Tiledistance>{Tiledistance{0,0,0,false}, Tiledistance{1,1,0,false},
                                                                          Tiledistance{2,2,0,false}, Tiledistance{3,3,0,false}}));
    // creates links in cubes (0,0,0,0), (0,0,1,0), (0,2,0,0), (0,2,1,0),
    //                        (0,-2,-2,-2), (0,-2,-1,-2), (0,0,-2,-2), (0,0,-1,-2)
    //                        (0,-4,-4,-4), (0,-4,-3,-4), (0,-2,-4,-4), (0,-2,-3,-4)
    std::vector<KmerOccurrence> occVector{
                                            KmerOccurrence{0,0,0,false,""}, KmerOccurrence{0,0,100,false,""}, KmerOccurrence{0,0,200,false,""},
                                            KmerOccurrence{1,1,0,false,""}, KmerOccurrence{1,1,110,false,""},
                                            KmerOccurrence{2,2,0,false,""}, KmerOccurrence{2,2,50,false,""},
                                            KmerOccurrence{3,3,0,false,""}
                                        };

    SECTION("With Hasse") {
        auto config = configBuilder.makeConfig();
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset(config, idMap);
        linkset.createRelevantLinks(occVector, 5, relevantCubes);
        tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> createdCubes;
//std::cout << "[DEBUG] -- relevant cubes " << relevantCubes << std::endl;
//std::cout << "[DEBUG] -- linkset " << linkset.linkset() << std::endl;
        for (auto&& elem : linkset.linkset()) {
            auto c = std::make_shared<Cube>(elem.first, config->tileSize());
            // only relevant links, possibly inducing lower-dim relevant cubes
            if (relevantCubes.find(c) == relevantCubes.end()) {
                Link l(std::vector<KmerOccurrence>{elem.first.occurrence(0), elem.first.occurrence(1)}, 5);
                auto c2 = std::make_shared<Cube>(l, config->tileSize());
                REQUIRE(relevantCubes.find(c2) != relevantCubes.end());
                createdCubes.insert(c2);
            } else {
                createdCubes.insert(c);
            }
        }
        REQUIRE(createdCubes == relevantCubes); // all relevant cubes covered
    }

    SECTION("No Hasse") {
        configBuilder.set(Hasse(false));
        auto config = configBuilder.makeConfig();
        Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset(config, idMap);
        auto fullCube = std::make_shared<Cube>(std::vector<Tiledistance>{Tiledistance{0,0,0,false}, Tiledistance{1,1,0,false},
                                                                         Tiledistance{2,2,0,false}, Tiledistance{3,3,0,false}});
        // links are created according to relevantCubes, processing function should not return incomplete occurrenceMap
        //   (i.e. where any genome has no occ), but if there is a relevantCube that is lower-dim and occurrences match,
        //   the according links will be created
        relevantCubes.clear();
        relevantCubes.insert(fullCube);
        linkset.createRelevantLinks(occVector, 5, relevantCubes);
        for (auto&& elem : linkset.linkset()) {
            Cube c(elem.first, config->tileSize());
            REQUIRE(c == *fullCube);
        }
    }
}
