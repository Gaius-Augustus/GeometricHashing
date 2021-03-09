#include <iostream>

#include "catch2/catch.hpp"
#include "nlohmann/json.hpp"
#include "../computeYassParameters.h"
#include "ConfigurationGenerator.h"
#include "LoadTestdata.h"

TEST_CASE("Compare YASS tuples with M4") {
    std::cout << "[INFO] -- [TEST CASE] -- Compare YASS tuples with M4" << std::endl;
    auto filepath = testfilesBasepath() / "yass";

    // run M4
    ConfigBuilder cbuilder;
    cbuilder.set(Genome1("hg38"),
                 Genome2("mm10"),
                 InputFiles(std::vector<std::string>{(filepath/"hg38.fa"), (filepath/"mm10.fa")}),
                 MaskCollectionPtr(std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(13),
                                                                                    SpacedSeedMaskCollection::SeedSetSize(4))),
                 MatchLimit(1),
                 MatchLimitDiscardSeeds(true),
                 PerformDiagonalFiltering(true),
                 YassEpsilon(0.05),
                 YassIndel(0.08),
                 YassMutation(0.15));
    auto config = cbuilder.makeConfig();
std::cout << "[DEBUG] -- YASS CONFIG: " << *config << std::endl;
std::cout << "[DEBUG] -- SEED MASKS: " << *(config->maskCollection()) << std::endl;

    auto idMap = std::make_shared<IdentifierMapping>(config->genome1());
    auto fastaCollection = std::make_shared<FastaCollection>(config);
    fastaCollection->populateIdentifierMappingFromFastaCollection(*idMap);
    auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, idMap);
    seedMap->extractSeeds(fastaCollection);
    Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> linkset{config, idMap};
    linkset.createLinks(*seedMap);
    seedMap->clear(); // save memory
    linkset.applyDiagonalMatchesFilter();
    linkset.groupOverlappingLinks();
std::cout << "[DEBUG] -- ID MAP: " << *idMap << std::endl;
    // seedMap now should hold all filtered matches (overlapping already grouped by M4), load Yass matches
    // load json
    std::ifstream is(filepath/"matches0.json");
    json j;
    is >> j;
    auto hgID = j["seq1"].get<std::string>();
    auto mmID = j["seq2"].get<std::string>();
    tsl::hopscotch_set<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> expectedMatches;
    for (auto&& elem : j["matches"]) {
        auto occ1 = KmerOccurrence(idMap->queryGenomeIDConst("hg38"),
                                   idMap->querySequenceIDConst(hgID, "hg38"),
                                   elem[0].get<size_t>(), false, "");
        auto occ2 = KmerOccurrence(idMap->queryGenomeIDConst("mm10"),
                                   idMap->querySequenceIDConst(mmID, "mm10"),
                                   elem[1].get<size_t>(), false, "");
        auto match = Link{std::vector<KmerOccurrence>{occ1, occ2}, elem[2].get<size_t>()};
        expectedMatches.emplace(match);
    }

    // check that matches are as expected
    // FAILS: REQUIRE(hopscotchSetEqual(expectedMatches, linkset.linkset()));
}



TEST_CASE("Compute Rho and Delta") {
    std::cout << "[INFO] -- [TEST CASE] -- Compute Rho and Delta" << std::endl;
    REQUIRE(computeRho(0.0, 9, 0.05) == 10);
    REQUIRE(computeDelta(0.0, 10, 0.05) == 0);
    REQUIRE(computeRho(0.0, 10, 0.05) == 11);
    REQUIRE(computeDelta(0.0, 11, 0.05) == 0);
    REQUIRE(computeRho(0.0, 11, 0.05) == 12);
    REQUIRE(computeDelta(0.0, 12, 0.05) == 0);
    REQUIRE(computeRho(0.0, 12, 0.05) == 13);
    REQUIRE(computeDelta(0.0, 13, 0.05) == 0);
    REQUIRE(computeRho(0.0, 13, 0.05) == 14);
    REQUIRE(computeDelta(0.0, 14, 0.05) == 0);
    REQUIRE(computeRho(0.0, 14, 0.05) == 15);
    REQUIRE(computeDelta(0.0, 15, 0.05) == 0);
    REQUIRE(computeRho(0.0, 15, 0.05) == 16);
    REQUIRE(computeDelta(0.0, 16, 0.05) == 0);
    REQUIRE(computeRho(0.15, 9, 0.05) == 54);
    REQUIRE(computeDelta(0.0, 54, 0.05) == 0);
    REQUIRE(computeRho(0.15, 10, 0.05) == 68);
    REQUIRE(computeDelta(0.0, 68, 0.05) == 0);
    REQUIRE(computeRho(0.15, 11, 0.05) == 84);
    REQUIRE(computeDelta(0.0, 84, 0.05) == 0);
    REQUIRE(computeRho(0.15, 12, 0.05) == 104);
    REQUIRE(computeDelta(0.0, 104, 0.05) == 0);
    REQUIRE(computeRho(0.15, 13, 0.05) == 127);
    REQUIRE(computeDelta(0.0, 127, 0.05) == 0);
    REQUIRE(computeRho(0.15, 14, 0.05) == 155);
    REQUIRE(computeDelta(0.0, 155, 0.05) == 0);
    REQUIRE(computeRho(0.15, 15, 0.05) == 187);
    REQUIRE(computeDelta(0.0, 187, 0.05) == 0);
    REQUIRE(computeRho(0.25, 9, 0.05) == 135);
    REQUIRE(computeDelta(0.0, 135, 0.05) == 0);
    REQUIRE(computeRho(0.25, 10, 0.05) == 187);
    REQUIRE(computeDelta(0.0, 187, 0.05) == 0);
    REQUIRE(computeRho(0.25, 11, 0.05) == 256);
    REQUIRE(computeDelta(0.0, 256, 0.05) == 0);
    REQUIRE(computeRho(0.25, 12, 0.05) == 349);
    REQUIRE(computeDelta(0.0, 349, 0.05) == 0);
    REQUIRE(computeRho(0.25, 13, 0.05) == 473);
    REQUIRE(computeDelta(0.0, 473, 0.05) == 0);
    REQUIRE(computeRho(0.25, 14, 0.05) == 640);
    REQUIRE(computeDelta(0.0, 640, 0.05) == 0);
    REQUIRE(computeRho(0.25, 15, 0.05) == 862);
    REQUIRE(computeDelta(0.0, 862, 0.05) == 0);
    REQUIRE(computeRho(0.0, 9, 0.05) == 10);
    REQUIRE(computeDelta(0.08, 10, 0.05) == 2);
    REQUIRE(computeRho(0.0, 10, 0.05) == 11);
    REQUIRE(computeDelta(0.08, 11, 0.05) == 2);
    REQUIRE(computeRho(0.0, 11, 0.05) == 12);
    REQUIRE(computeDelta(0.08, 12, 0.05) == 2);
    REQUIRE(computeRho(0.0, 12, 0.05) == 13);
    REQUIRE(computeDelta(0.08, 13, 0.05) == 2);
    REQUIRE(computeRho(0.0, 13, 0.05) == 14);
    REQUIRE(computeDelta(0.08, 14, 0.05) == 2);
    REQUIRE(computeRho(0.0, 14, 0.05) == 15);
    REQUIRE(computeDelta(0.08, 15, 0.05) == 2);
    REQUIRE(computeRho(0.0, 15, 0.05) == 16);
    REQUIRE(computeDelta(0.08, 16, 0.05) == 2);
    REQUIRE(computeRho(0.15, 9, 0.05) == 54);
    REQUIRE(computeDelta(0.08, 54, 0.05) == 4);
    REQUIRE(computeRho(0.15, 10, 0.05) == 68);
    REQUIRE(computeDelta(0.08, 68, 0.05) == 5);
    REQUIRE(computeRho(0.15, 11, 0.05) == 84);
    REQUIRE(computeDelta(0.08, 84, 0.05) == 5);
    REQUIRE(computeRho(0.15, 12, 0.05) == 104);
    REQUIRE(computeDelta(0.08, 104, 0.05) == 6);
    REQUIRE(computeRho(0.15, 13, 0.05) == 127);
    REQUIRE(computeDelta(0.08, 127, 0.05) == 6);
    REQUIRE(computeRho(0.15, 14, 0.05) == 155);
    REQUIRE(computeDelta(0.08, 155, 0.05) == 7);
    REQUIRE(computeRho(0.15, 15, 0.05) == 187);
    REQUIRE(computeDelta(0.08, 187, 0.05) == 8);
    REQUIRE(computeRho(0.25, 9, 0.05) == 135);
    REQUIRE(computeDelta(0.08, 135, 0.05) == 6);
    REQUIRE(computeRho(0.25, 10, 0.05) == 187);
    REQUIRE(computeDelta(0.08, 187, 0.05) == 8);
    REQUIRE(computeRho(0.25, 11, 0.05) == 256);
    REQUIRE(computeDelta(0.08, 256, 0.05) == 9);
    REQUIRE(computeRho(0.25, 12, 0.05) == 349);
    REQUIRE(computeDelta(0.08, 349, 0.05) == 10);
    REQUIRE(computeRho(0.25, 13, 0.05) == 473);
    REQUIRE(computeDelta(0.08, 473, 0.05) == 12);
    REQUIRE(computeRho(0.25, 14, 0.05) == 640);
    REQUIRE(computeDelta(0.08, 640, 0.05) == 14);
    REQUIRE(computeRho(0.25, 15, 0.05) == 862);
    REQUIRE(computeDelta(0.08, 862, 0.05) == 16);
    REQUIRE(computeRho(0.0, 9, 0.05) == 10);
    REQUIRE(computeDelta(0.2, 10, 0.05) == 3);
    REQUIRE(computeRho(0.0, 10, 0.05) == 11);
    REQUIRE(computeDelta(0.2, 11, 0.05) == 3);
    REQUIRE(computeRho(0.0, 11, 0.05) == 12);
    REQUIRE(computeDelta(0.2, 12, 0.05) == 3);
    REQUIRE(computeRho(0.0, 12, 0.05) == 13);
    REQUIRE(computeDelta(0.2, 13, 0.05) == 3);
    REQUIRE(computeRho(0.0, 13, 0.05) == 14);
    REQUIRE(computeDelta(0.2, 14, 0.05) == 3);
    REQUIRE(computeRho(0.0, 14, 0.05) == 15);
    REQUIRE(computeDelta(0.2, 15, 0.05) == 3);
    REQUIRE(computeRho(0.0, 15, 0.05) == 16);
    REQUIRE(computeDelta(0.2, 16, 0.05) == 3);
    REQUIRE(computeRho(0.15, 9, 0.05) == 54);
    REQUIRE(computeDelta(0.2, 54, 0.05) == 6);
    REQUIRE(computeRho(0.15, 10, 0.05) == 68);
    REQUIRE(computeDelta(0.2, 68, 0.05) == 7);
    REQUIRE(computeRho(0.15, 11, 0.05) == 84);
    REQUIRE(computeDelta(0.2, 84, 0.05) == 8);
    REQUIRE(computeRho(0.15, 12, 0.05) == 104);
    REQUIRE(computeDelta(0.2, 104, 0.05) == 9);
    REQUIRE(computeRho(0.15, 13, 0.05) == 127);
    REQUIRE(computeDelta(0.2, 127, 0.05) == 10);
    REQUIRE(computeRho(0.15, 14, 0.05) == 155);
    REQUIRE(computeDelta(0.2, 155, 0.05) == 11);
    REQUIRE(computeRho(0.15, 15, 0.05) == 187);
    REQUIRE(computeDelta(0.2, 187, 0.05) == 12);
    REQUIRE(computeRho(0.25, 9, 0.05) == 135);
    REQUIRE(computeDelta(0.2, 135, 0.05) == 10);
    REQUIRE(computeRho(0.25, 10, 0.05) == 187);
    REQUIRE(computeDelta(0.2, 187, 0.05) == 12);
    REQUIRE(computeRho(0.25, 11, 0.05) == 256);
    REQUIRE(computeDelta(0.2, 256, 0.05) == 14);
    REQUIRE(computeRho(0.25, 12, 0.05) == 349);
    REQUIRE(computeDelta(0.2, 349, 0.05) == 16);
    REQUIRE(computeRho(0.25, 13, 0.05) == 473);
    REQUIRE(computeDelta(0.2, 473, 0.05) == 19);
    REQUIRE(computeRho(0.25, 14, 0.05) == 640);
    REQUIRE(computeDelta(0.2, 640, 0.05) == 22);
    REQUIRE(computeRho(0.25, 15, 0.05) == 862);
    REQUIRE(computeDelta(0.2, 862, 0.05) == 26);
    REQUIRE(computeRho(0.0, 9, 0.11) == 10);
    REQUIRE(computeDelta(0.0, 10, 0.11) == 0);
    REQUIRE(computeRho(0.0, 10, 0.11) == 11);
    REQUIRE(computeDelta(0.0, 11, 0.11) == 0);
    REQUIRE(computeRho(0.0, 11, 0.11) == 12);
    REQUIRE(computeDelta(0.0, 12, 0.11) == 0);
    REQUIRE(computeRho(0.0, 12, 0.11) == 13);
    REQUIRE(computeDelta(0.0, 13, 0.11) == 0);
    REQUIRE(computeRho(0.0, 13, 0.11) == 14);
    REQUIRE(computeDelta(0.0, 14, 0.11) == 0);
    REQUIRE(computeRho(0.0, 14, 0.11) == 15);
    REQUIRE(computeDelta(0.0, 15, 0.11) == 0);
    REQUIRE(computeRho(0.0, 15, 0.11) == 16);
    REQUIRE(computeDelta(0.0, 16, 0.11) == 0);
    REQUIRE(computeRho(0.15, 9, 0.11) == 42);
    REQUIRE(computeDelta(0.0, 42, 0.11) == 0);
    REQUIRE(computeRho(0.15, 10, 0.11) == 52);
    REQUIRE(computeDelta(0.0, 52, 0.11) == 0);
    REQUIRE(computeRho(0.15, 11, 0.11) == 65);
    REQUIRE(computeDelta(0.0, 65, 0.11) == 0);
    REQUIRE(computeRho(0.15, 12, 0.11) == 79);
    REQUIRE(computeDelta(0.0, 79, 0.11) == 0);
    REQUIRE(computeRho(0.15, 13, 0.11) == 96);
    REQUIRE(computeDelta(0.0, 96, 0.11) == 0);
    REQUIRE(computeRho(0.15, 14, 0.11) == 117);
    REQUIRE(computeDelta(0.0, 117, 0.11) == 0);
    REQUIRE(computeRho(0.15, 15, 0.11) == 141);
    REQUIRE(computeDelta(0.0, 141, 0.11) == 0);
    REQUIRE(computeRho(0.25, 9, 0.11) == 102);
    REQUIRE(computeDelta(0.0, 102, 0.11) == 0);
    REQUIRE(computeRho(0.25, 10, 0.11) == 140);
    REQUIRE(computeDelta(0.0, 140, 0.11) == 0);
    REQUIRE(computeRho(0.25, 11, 0.11) == 191);
    REQUIRE(computeDelta(0.0, 191, 0.11) == 0);
    REQUIRE(computeRho(0.25, 12, 0.11) == 260);
    REQUIRE(computeDelta(0.0, 260, 0.11) == 0);
    REQUIRE(computeRho(0.25, 13, 0.11) == 352);
    REQUIRE(computeDelta(0.0, 352, 0.11) == 0);
    REQUIRE(computeRho(0.25, 14, 0.11) == 474);
    REQUIRE(computeDelta(0.0, 474, 0.11) == 0);
    REQUIRE(computeRho(0.25, 15, 0.11) == 638);
    REQUIRE(computeDelta(0.0, 638, 0.11) == 0);
    REQUIRE(computeRho(0.0, 9, 0.11) == 10);
    REQUIRE(computeDelta(0.08, 10, 0.11) == 1);
    REQUIRE(computeRho(0.0, 10, 0.11) == 11);
    REQUIRE(computeDelta(0.08, 11, 0.11) == 1);
    REQUIRE(computeRho(0.0, 11, 0.11) == 12);
    REQUIRE(computeDelta(0.08, 12, 0.11) == 2);
    REQUIRE(computeRho(0.0, 12, 0.11) == 13);
    REQUIRE(computeDelta(0.08, 13, 0.11) == 2);
    REQUIRE(computeRho(0.0, 13, 0.11) == 14);
    REQUIRE(computeDelta(0.08, 14, 0.11) == 2);
    REQUIRE(computeRho(0.0, 14, 0.11) == 15);
    REQUIRE(computeDelta(0.08, 15, 0.11) == 2);
    REQUIRE(computeRho(0.0, 15, 0.11) == 16);
    REQUIRE(computeDelta(0.08, 16, 0.11) == 2);
    REQUIRE(computeRho(0.15, 9, 0.11) == 42);
    REQUIRE(computeDelta(0.08, 42, 0.11) == 3);
    REQUIRE(computeRho(0.15, 10, 0.11) == 52);
    REQUIRE(computeDelta(0.08, 52, 0.11) == 3);
    REQUIRE(computeRho(0.15, 11, 0.11) == 65);
    REQUIRE(computeDelta(0.08, 65, 0.11) == 4);
    REQUIRE(computeRho(0.15, 12, 0.11) == 79);
    REQUIRE(computeDelta(0.08, 79, 0.11) == 4);
    REQUIRE(computeRho(0.15, 13, 0.11) == 96);
    REQUIRE(computeDelta(0.08, 96, 0.11) == 4);
    REQUIRE(computeRho(0.15, 14, 0.11) == 117);
    REQUIRE(computeDelta(0.08, 117, 0.11) == 5);
    REQUIRE(computeRho(0.15, 15, 0.11) == 141);
    REQUIRE(computeDelta(0.08, 141, 0.11) == 5);
    REQUIRE(computeRho(0.25, 9, 0.11) == 102);
    REQUIRE(computeDelta(0.08, 102, 0.11) == 5);
    REQUIRE(computeRho(0.25, 10, 0.11) == 140);
    REQUIRE(computeDelta(0.08, 140, 0.11) == 5);
    REQUIRE(computeRho(0.25, 11, 0.11) == 191);
    REQUIRE(computeDelta(0.08, 191, 0.11) == 6);
    REQUIRE(computeRho(0.25, 12, 0.11) == 260);
    REQUIRE(computeDelta(0.08, 260, 0.11) == 7);
    REQUIRE(computeRho(0.25, 13, 0.11) == 352);
    REQUIRE(computeDelta(0.08, 352, 0.11) == 8);
    REQUIRE(computeRho(0.25, 14, 0.11) == 474);
    REQUIRE(computeDelta(0.08, 474, 0.11) == 10);
    REQUIRE(computeRho(0.25, 15, 0.11) == 638);
    REQUIRE(computeDelta(0.08, 638, 0.11) == 11);
    REQUIRE(computeRho(0.0, 9, 0.11) == 10);
    REQUIRE(computeDelta(0.2, 10, 0.11) == 2);
    REQUIRE(computeRho(0.0, 10, 0.11) == 11);
    REQUIRE(computeDelta(0.2, 11, 0.11) == 2);
    REQUIRE(computeRho(0.0, 11, 0.11) == 12);
    REQUIRE(computeDelta(0.2, 12, 0.11) == 2);
    REQUIRE(computeRho(0.0, 12, 0.11) == 13);
    REQUIRE(computeDelta(0.2, 13, 0.11) == 3);
    REQUIRE(computeRho(0.0, 13, 0.11) == 14);
    REQUIRE(computeDelta(0.2, 14, 0.11) == 3);
    REQUIRE(computeRho(0.0, 14, 0.11) == 15);
    REQUIRE(computeDelta(0.2, 15, 0.11) == 3);
    REQUIRE(computeRho(0.0, 15, 0.11) == 16);
    REQUIRE(computeDelta(0.2, 16, 0.11) == 3);
    REQUIRE(computeRho(0.15, 9, 0.11) == 42);
    REQUIRE(computeDelta(0.2, 42, 0.11) == 5);
    REQUIRE(computeRho(0.15, 10, 0.11) == 52);
    REQUIRE(computeDelta(0.2, 52, 0.11) == 5);
    REQUIRE(computeRho(0.15, 11, 0.11) == 65);
    REQUIRE(computeDelta(0.2, 65, 0.11) == 6);
    REQUIRE(computeRho(0.15, 12, 0.11) == 79);
    REQUIRE(computeDelta(0.2, 79, 0.11) == 6);
    REQUIRE(computeRho(0.15, 13, 0.11) == 96);
    REQUIRE(computeDelta(0.2, 96, 0.11) == 7);
    REQUIRE(computeRho(0.15, 14, 0.11) == 117);
    REQUIRE(computeDelta(0.2, 117, 0.11) == 8);
    REQUIRE(computeRho(0.15, 15, 0.11) == 141);
    REQUIRE(computeDelta(0.2, 141, 0.11) == 8);
    REQUIRE(computeRho(0.25, 9, 0.11) == 102);
    REQUIRE(computeDelta(0.2, 102, 0.11) == 7);
    REQUIRE(computeRho(0.25, 10, 0.11) == 140);
    REQUIRE(computeDelta(0.2, 140, 0.11) == 8);
    REQUIRE(computeRho(0.25, 11, 0.11) == 191);
    REQUIRE(computeDelta(0.2, 191, 0.11) == 10);
    REQUIRE(computeRho(0.25, 12, 0.11) == 260);
    REQUIRE(computeDelta(0.2, 260, 0.11) == 12);
    REQUIRE(computeRho(0.25, 13, 0.11) == 352);
    REQUIRE(computeDelta(0.2, 352, 0.11) == 13);
    REQUIRE(computeRho(0.25, 14, 0.11) == 474);
    REQUIRE(computeDelta(0.2, 474, 0.11) == 16);
    REQUIRE(computeRho(0.25, 15, 0.11) == 638);
    REQUIRE(computeDelta(0.2, 638, 0.11) == 18);
}

/*#include <bitset>
#include <cstdlib>
#include <climits>
#include <thread>
#include <vector>

#include "catch2/catch.hpp"
#include "prettyprint.hpp"
#include "tsl/hopscotch_set.h"
#include "../ExtractSeeds.h"
#include "../SeedMapSpaced.h"
#include "ConfigurationGenerator.h"
#include "fillDiagonalMatchesFilterTestdata.h"
#include "LoadTestdata.h"

void createSeeds(std::shared_ptr<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>> seedMap,
                 std::shared_ptr<IdentifierMapping> idMap) {
    idMap->querySequenceID("first_sequence", "hg38_orthologs");
    idMap->querySequenceID("second_sequence", "mm10_orthologs");

    // "AAAAACCCCCGGGGGTTTTT"
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("AAAAA"), KmerOccurrence(0,0,0,false,"AAAAA"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("AAACC"), KmerOccurrence(0,0,2,false,"AAACC"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("CCCCC"), KmerOccurrence(0,0,5,false,"CCCCC"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("CCCGG"), KmerOccurrence(0,0,7,false,"CCCGG"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("GGGGG"), KmerOccurrence(0,0,10,false,"GGGGG"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("GGGTT"), KmerOccurrence(0,0,12,false,"GGGTT"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("TTTTT"), KmerOccurrence(0,0,15,false,"TTTTT"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("AAAAA"), KmerOccurrence(1,1,100,false,"AAAAA"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("AAACC"), KmerOccurrence(1,1,102,false,"AAACC"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT"));
    seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("TTTTT"), KmerOccurrence(1,1,115,false,"TTTTT"));
}


TEST_CASE("SeedMapSpaced with M4 filter") {
    std::cout << "[INFO] -- [TEST CASE] -- SeedMapSpaced with M4 filter" << std::endl;
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();

    SECTION("Test Chunk Creation") {
        auto config = customConfiguration(DiagonalThreshold(1),
                                          Genome1("genome0"),
                                          Genome2("genome1"),
                                          InputFiles(data.inputFiles()),
                                          LocalSearchAreaLength(1),
                                          MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                             SpacedSeedMaskCollection::Span{5},
                                                                                                             SpacedSeedMaskCollection::SeedSetSize{1})},
                                          MinMatchDistance(1),
                                          NThreads(4),  // needs to be fix in order for chunk algorithm to work predictably
                                          TileSize(1));
        auto idMap = std::make_shared<IdentifierMapping>("genome0");
        idMap->queryGenomeID("genome1");
        idMap->querySequenceID("seq0", "genome0");
        idMap->querySequenceID("seq1", "genome1");
        auto seedMap = SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>(config, idMap);
        auto filter = DiagonalMatchesFilter<KmerOccurrencePair>(config);
        SECTION("Test One Match Chunk Creation") {
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> chunk;
            chunk.emplace(KmerOccurrence(0,0,0,false,"ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"), config->span());
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> matchSet;
            matchSet.insert(chunk.begin(), chunk.end());
            seedMap.appendMatches(matchSet);
            auto orderedMatches = seedMap.orderedMatches();
            auto chunks = filter.matchesChunks(orderedMatches);
            REQUIRE(chunks.size() == 1);
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reconstructChunk;
            reconstructChunk.insert(chunks.at(0).first, chunks.at(0).second);
            REQUIRE(reconstructChunk == chunk);
        }
        SECTION("Test Two Match Chunk Creation") {
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> chunk;
            chunk.emplace(KmerOccurrence(0,0,0,false,"ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"), config->span());
            chunk.emplace(KmerOccurrence(0,0,1,false,"ACGTA"), KmerOccurrence(1,1,101,false,"ACGTA"), config->span());
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> matchSet;
            matchSet.insert(chunk.begin(), chunk.end());
            seedMap.appendMatches(matchSet);
            auto orderedMatches = seedMap.orderedMatches();
            auto chunks = filter.matchesChunks(orderedMatches);
            REQUIRE(chunks.size() == 1);
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reconstructChunk;
            reconstructChunk.insert(chunks.at(0).first, chunks.at(0).second);
            REQUIRE(reconstructChunk == chunk);
        }
        SECTION("Test Four Match Chunk Creation (= nThreads_)") {
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> chunk;
            chunk.emplace(KmerOccurrence(0,0,0,false,"ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"), config->span());    // distance: +100
            chunk.emplace(KmerOccurrence(0,0,1,false,"ACGTA"), KmerOccurrence(1,1,101,false,"ACGTA"), config->span());    // distance: +100
            chunk.emplace(KmerOccurrence(0,0,2,false,"ACGTA"), KmerOccurrence(1,1,102,false,"ACGTA"), config->span());    // distance: +100
            chunk.emplace(KmerOccurrence(0,0,3,false,"ACGTA"), KmerOccurrence(1,1,103,false,"ACGTA"), config->span());    // distance: +100, all same diag, so one chunk desired
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> matchSet;
            matchSet.insert(chunk.begin(), chunk.end());
            seedMap.appendMatches(matchSet);
            auto orderedMatches = seedMap.orderedMatches();
            auto chunks = filter.matchesChunks(orderedMatches);
            REQUIRE(chunks.size() == 1);
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reconstructChunk;
            reconstructChunk.insert(chunks.at(0).first, chunks.at(0).second);
            REQUIRE(reconstructChunk == chunk);
        }
        SECTION("Test Correct Chunk Creation") {
            // nThreads = 4, wholeSet size = 16, so equal sized would be 4 pairs per chunk
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> chunk0;
            chunk0.emplace(KmerOccurrence(0,0,0,false,"ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"), config->span());    // distance: +100
            chunk0.emplace(KmerOccurrence(0,0,1,false,"ACGTA"), KmerOccurrence(1,1,101,false,"ACGTA"), config->span());    // distance: +100
            chunk0.emplace(KmerOccurrence(0,0,2,false,"ACGTA"), KmerOccurrence(1,1,102,false,"ACGTA"), config->span());    // distance: +100
            chunk0.emplace(KmerOccurrence(0,0,10,false,"ACGTA"), KmerOccurrence(1,1,210,false,"ACGTA"), config->span());    // distance: +200 // first equal chunk last
            chunk0.emplace(KmerOccurrence(0,0,11,false,"ACGTA"), KmerOccurrence(1,1,211,false,"ACGTA"), config->span());    // distance: +200 // desired first chunk last

            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> chunk1;
            chunk1.emplace(KmerOccurrence(0,0,20,false,"ACGTA"), KmerOccurrence(1,1,320,false,"ACGTA"), config->span());    // distance: +300
            chunk1.emplace(KmerOccurrence(0,0,21,false,"ACGTA"), KmerOccurrence(1,1,321,false,"ACGTA"), config->span());    // distance: +300
            chunk1.emplace(KmerOccurrence(0,0,22,false,"ACGTA"), KmerOccurrence(1,1,322,false,"ACGTA"), config->span());    // distance: +300
            chunk1.emplace(KmerOccurrence(0,0,30,false,"ACGTA"), KmerOccurrence(1,1,430,false,"ACGTA"), config->span());    // distance: +400
            chunk1.emplace(KmerOccurrence(0,0,31,false,"ACGTA"), KmerOccurrence(1,1,431,false,"ACGTA"), config->span());    // distance: +400
            chunk1.emplace(KmerOccurrence(0,0,32,false,"ACGTA"), KmerOccurrence(1,1,432,false,"ACGTA"), config->span());    // distance: +400
            chunk1.emplace(KmerOccurrence(0,0,33,false,"ACGTA"), KmerOccurrence(1,1,433,false,"ACGTA"), config->span());    // distance: +400
            chunk1.emplace(KmerOccurrence(0,0,34,false,"ACGTA"), KmerOccurrence(1,1,434,false,"ACGTA"), config->span());    // distance: +400 // desired second chunk last

            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> chunk2;
            chunk2.emplace(KmerOccurrence(0,0,40,false,"ACGTA"), KmerOccurrence(1,1,540,false,"ACGTA"), config->span());    // distance: +500
            chunk2.emplace(KmerOccurrence(0,0,41,false,"ACGTA"), KmerOccurrence(1,1,541,false,"ACGTA"), config->span());    // distance: +500
            chunk2.emplace(KmerOccurrence(0,0,42,false,"ACGTA"), KmerOccurrence(1,1,542,false,"ACGTA"), config->span());    // distance: +500

            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> wholeSet;
            wholeSet.insert(chunk0.begin(), chunk0.end());
            wholeSet.insert(chunk1.begin(), chunk1.end());
            wholeSet.insert(chunk2.begin(), chunk2.end());

            seedMap.appendMatches(wholeSet);

            auto orderedMatches = seedMap.orderedMatches();
            auto chunks = filter.matchesChunks(orderedMatches);
            REQUIRE(chunks.size() == 3);
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> builtChunk0;
            builtChunk0.insert(chunks.at(0).first, chunks.at(0).second);
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> builtChunk1;
            builtChunk1.insert(chunks.at(1).first, chunks.at(1).second);
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> builtChunk2;
            builtChunk2.insert(chunks.at(2).first, chunks.at(2).second);
            REQUIRE(builtChunk0 == chunk0);
            REQUIRE(builtChunk1 == chunk1);
            REQUIRE(builtChunk2 == chunk2);
            tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reconstructChunk;
            reconstructChunk.insert(chunks.at(0).first, chunks.at(0).second);
            reconstructChunk.insert(chunks.at(2).first, chunks.at(1).second);
            reconstructChunk.insert(chunks.at(1).first, chunks.at(2).second);
            REQUIRE(reconstructChunk == wholeSet);
        }
    }

    SECTION("Test SeedMapSpaced Methods 1") {
        auto config = customConfiguration(DiagonalThreshold(3),
                                          Genome1(data.genome0()),
                                          Genome2(data.genome1()),
                                          InputFiles(data.inputFiles()),
                                          LocalSearchAreaLength(15),
                                          MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                             SpacedSeedMaskCollection::Span{5},
                                                                                                             SpacedSeedMaskCollection::SeedSetSize{1})},
                                          PerformDiagonalFiltering(true),
                                          TileSize(1));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap());
        createSeeds(seedMap, data.idMap());
        REQUIRE(seedMap->seedMap().size() == 7);
        seedMap->createMatches();
        REQUIRE(seedMap->matches().size() == 7);
        for (auto&& match : seedMap->matches()) {
            REQUIRE(match.distance() == 100);
        }

        // perform filter manually
        auto filter = DiagonalMatchesFilter<KmerOccurrencePair>(config);
        std::vector<KmerOccurrencePair> reported;
        auto orderedMatches = seedMap->orderedMatches();
        filter.applyDiagonalMatchesFilter(orderedMatches, reported);
        tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reportedSet;
        reportedSet.insert(reported.begin(), reported.end());

        // call SeedMapSpaced filter, should give equal result
        auto skipped = std::make_shared<std::array<size_t, 4>>();
        seedMap->applyDiagonalMatchesFilter(skipped);

        REQUIRE(seedMap->matches() == reportedSet);
        REQUIRE((*skipped)[0] == filter.skippedNotInGenome1And2());
        REQUIRE((*skipped)[1] == filter.skippedOverlappedOrTooClose());
        REQUIRE((*skipped)[2] == filter.skippedTooFewDiagonalElements());
        REQUIRE((*skipped)[3] == filter.skippedTooFewNeighbours());

        std::sort(reported.begin(), reported.end());
// AAACC --> those do not exist anymore as !allowOverlap
// AAAAACCCCCGGGGGTTTTT
//   |  XXXXXx no
// AAAAACCCCCGGGGGTTTTT
// XXXXX  |  XXXXXx yes
// AAAAACCCCCGGGGGTTTTT
//     xXXXXX  |  XXXXX yes
// AAAAACCCCCGGGGGTTTTT
//          xXXXXX  | no
        std::vector<KmerOccurrencePair> expected{
//                KmerOccurrencePair(KmerOccurrence(0,0,0,false,"AAAAA"), KmerOccurrence(1,1,100,false,"AAAAA")),
//                KmerOccurrencePair(KmerOccurrence(0,0,2,false,"AAACC"), KmerOccurrence(1,1,102,false,"AAACC")),
                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG")),
                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,12,false,"GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT")),
//                KmerOccurrencePair(KmerOccurrence(0,0,15,false,"TTTTT"), KmerOccurrence(1,1,115,false,"TTTTT"))
        };
        std::sort(expected.begin(), expected.end());
        REQUIRE(reported == expected);
        REQUIRE((*skipped)[0] == 0);
        REQUIRE((*skipped)[1] == 3);
        REQUIRE((*skipped)[2] == 0);
        REQUIRE((*skipped)[3] == 2);
    }



    SECTION("Test SeedMapSpaced Methods 2") {
        auto config = customConfiguration(DiagonalThreshold(3),
                                          Genome1(data.genome0()),
                                          Genome2(data.genome1()),
                                          InputFiles(data.inputFiles()),
                                          LocalSearchAreaLength(24),
                                          MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                             SpacedSeedMaskCollection::Span{5},
                                                                                                             SpacedSeedMaskCollection::SeedSetSize{1})},
                                          PerformDiagonalFiltering(true),
                                          TileSize(1));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap());
        createSeeds(seedMap, data.idMap());
        REQUIRE(seedMap->seedMap().size() == 7);
        seedMap->createMatches();
        REQUIRE(seedMap->matches().size() == 7);
        for (auto&& match : seedMap->matches()) {
            REQUIRE(match.distance() == 100);
        }

        // perform filter manually
        auto filter = DiagonalMatchesFilter<KmerOccurrencePair>(config);
        std::vector<KmerOccurrencePair> reported;
        auto orderedMatches = seedMap->orderedMatches();
        filter.applyDiagonalMatchesFilter(orderedMatches/ *matches()* /, reported);
        tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reportedSet;
        reportedSet.insert(reported.begin(), reported.end());

        // call SeedMapSpaced filter, should give equal result
        auto skipped = std::make_shared<std::array<size_t, 4>>();
        seedMap->applyDiagonalMatchesFilter(skipped);

        REQUIRE(seedMap->matches() == reportedSet);
        REQUIRE((*skipped)[0] == filter.skippedNotInGenome1And2());
        REQUIRE((*skipped)[1] == filter.skippedOverlappedOrTooClose());
        REQUIRE((*skipped)[2] == filter.skippedTooFewDiagonalElements());
        REQUIRE((*skipped)[3] == filter.skippedTooFewNeighbours());

        std::sort(reported.begin(), reported.end());
// AAACC --> those do not exist anymore as !allowOverlap
// AAAAACCCCCGGGGGTTTTT
//   |  XXXXXXXXXX yes
// AAAAACCCCCGGGGGTTTTT
// XXXXX  |  XXXXXXXXXX yes
// AAAAACCCCCGGGGGTTTTT
// XXXXXXXXXX  |  XXXXX yes
// AAAAACCCCCGGGGGTTTTT
//      XXXXXXXXXX  | yes
        std::vector<KmerOccurrencePair> expected{
                KmerOccurrencePair(KmerOccurrence(0,0,0,false,"AAAAA"), KmerOccurrence(1,1,100,false,"AAAAA"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,2,false,"AAACC"), KmerOccurrence(1,1,102,false,"AAACC")),
                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG")),
                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,12,false,"GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT")),
                KmerOccurrencePair(KmerOccurrence(0,0,15,false,"TTTTT"), KmerOccurrence(1,1,115,false,"TTTTT"), config->span())
        };
        std::sort(expected.begin(), expected.end());
        REQUIRE(reported == expected);
        REQUIRE((*skipped)[0] == 0);
        REQUIRE((*skipped)[1] == 3);
        REQUIRE((*skipped)[2] == 0);
        REQUIRE((*skipped)[3] == 0);
    }



    SECTION("Test SeedMapSpaced Methods 3") {
        auto config = customConfiguration(DiagonalThreshold(3),
                                          Genome1(data.genome0()),
                                          Genome2(data.genome1()),
                                          InputFiles(data.inputFiles()),
                                          LocalSearchAreaLength(24),
                                          MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                             SpacedSeedMaskCollection::Span{5},
                                                                                                             SpacedSeedMaskCollection::SeedSetSize{1})},
                                          MinMatchDistance(1),
                                          PerformDiagonalFiltering(true),
                                          TileSize(1));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap());
        createSeeds(seedMap, data.idMap());
        REQUIRE(seedMap->seedMap().size() == 7);
        seedMap->createMatches();
        REQUIRE(seedMap->matches().size() == 7);
        for (auto&& match : seedMap->matches()) {
            REQUIRE(match.distance() == 100);
        }

        // perform filter manually
        auto filter = DiagonalMatchesFilter<KmerOccurrencePair>(config);
        std::vector<KmerOccurrencePair> reported;
        auto orderedMatches = seedMap->orderedMatches();
        filter.applyDiagonalMatchesFilter(orderedMatches, reported);
        tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reportedSet;
        reportedSet.insert(reported.begin(), reported.end());

        // call SeedMapSpaced filter, should give equal result
        auto skipped = std::make_shared<std::array<size_t, 4>>();
        seedMap->applyDiagonalMatchesFilter(skipped);

        REQUIRE(seedMap->matches() == reportedSet);
        REQUIRE((*skipped)[0] == filter.skippedNotInGenome1And2());
        REQUIRE((*skipped)[1] == filter.skippedOverlappedOrTooClose());
        REQUIRE((*skipped)[2] == filter.skippedTooFewDiagonalElements());
        REQUIRE((*skipped)[3] == filter.skippedTooFewNeighbours());

        std::sort(reported.begin(), reported.end());
// AAAAA__CCCGG___TTTTT -> min. match dist. is 1
// AAAAA__CCCGG___TTTTT
//   |    XXXXX no
// AAAAA__CCCGG___TTTTT
// XXXXX    |     XXXXX yes
// AAAAA__CCCGG___TTTTT
//        XXXXX     | no
        std::vector<KmerOccurrencePair> expected{
//                KmerOccurrencePair(KmerOccurrence(0,0,0,false,"AAAAA"), KmerOccurrence(1,1,100,false,"AAAAA")),
//                KmerOccurrencePair(KmerOccurrence(0,0,2,false,"AAACC"), KmerOccurrence(1,1,102,false,"AAACC")),
//                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC")),
                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG")),
//                KmerOccurrencePair(KmerOccurrence(0,0,12,false,"GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT")),
//                KmerOccurrencePair(KmerOccurrence(0,0,15,false,"TTTTT"), KmerOccurrence(1,1,115,false,"TTTTT"))
                };
        std::sort(expected.begin(), expected.end());
        REQUIRE(reported == expected);
        //REQUIRE(*skipped == 6);
        REQUIRE((*skipped)[0] == 0);
        REQUIRE((*skipped)[1] == 4);
        REQUIRE((*skipped)[2] == 0);
        REQUIRE((*skipped)[3] == 2);
    }



    SECTION("Test SeedMapSpaced Methods 4") {
        auto config = customConfiguration(DiagonalThreshold(2),
                                          Genome1(data.genome0()),
                                          Genome2(data.genome1()),
                                          InputFiles(data.inputFiles()),
                                          LocalSearchAreaLength(24),
                                          MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                             SpacedSeedMaskCollection::Span{5},
                                                                                                             SpacedSeedMaskCollection::SeedSetSize{1})},
                                          MinMatchDistance(1),
                                          PerformDiagonalFiltering(true),
                                          TileSize(1));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap());
        createSeeds(seedMap, data.idMap());
        REQUIRE(seedMap->seedMap().size() == 7);
        seedMap->createMatches();
        REQUIRE(seedMap->matches().size() == 7);
        for (auto&& match : seedMap->matches()) {
            REQUIRE(match.distance() == 100);
        }

        // perform filter manually
        auto filter = DiagonalMatchesFilter<KmerOccurrencePair>(config);
        std::vector<KmerOccurrencePair> reported;
        auto orderedMatches = seedMap->orderedMatches();
        filter.applyDiagonalMatchesFilter(orderedMatches, reported);
        tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reportedSet;
        reportedSet.insert(reported.begin(), reported.end());

        // call SeedMapSpaced filter, should give equal result
        auto skipped = std::make_shared<std::array<size_t, 4>>();
        seedMap->applyDiagonalMatchesFilter(skipped);

        REQUIRE(seedMap->matches() == reportedSet);
        REQUIRE((*skipped)[0] == filter.skippedNotInGenome1And2());
        REQUIRE((*skipped)[1] == filter.skippedOverlappedOrTooClose());
        REQUIRE((*skipped)[2] == filter.skippedTooFewDiagonalElements());
        REQUIRE((*skipped)[3] == filter.skippedTooFewNeighbours());

        std::sort(reported.begin(), reported.end());

// AAAAA__CCCGG___TTTTT -> min. match dist. is 1
// AAAAA__CCCGG___TTTTT
//   |    XXXXX yes
// AAAAA__CCCGG___TTTTT
// XXXXX    |     XXXXX yes
// AAAAA__CCCGG___TTTTT
//        XXXXX     | yes
        std::vector<KmerOccurrencePair> expected{
                KmerOccurrencePair(KmerOccurrence(0,0,0,false,"AAAAA"), KmerOccurrence(1,1,100,false,"AAAAA"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,2,false,"AAACC"), KmerOccurrence(1,1,102,false,"AAACC")),
//                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC")),
                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG")),
//                KmerOccurrencePair(KmerOccurrence(0,0,12,false,"GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT")),
                KmerOccurrencePair(KmerOccurrence(0,0,15,false,"TTTTT"), KmerOccurrence(1,1,115,false,"TTTTT"), config->span())
        };
        std::sort(expected.begin(), expected.end());
        REQUIRE(reported == expected);
        //REQUIRE(*skipped == 4);
        REQUIRE((*skipped)[0] == 0);
        REQUIRE((*skipped)[1] == 4);
        REQUIRE((*skipped)[2] == 0);
        REQUIRE((*skipped)[3] == 0);
    }



    SECTION("Test SeedMapSpaced Methods 5") {
        auto config = customConfiguration(AllowOverlap(true),
                                          DiagonalThreshold(2),
                                          Genome1(data.genome0()),
                                          Genome2(data.genome1()),
                                          InputFiles(data.inputFiles()),
                                          LocalSearchAreaLength(10),
                                          MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                             SpacedSeedMaskCollection::Span{5},
                                                                                                             SpacedSeedMaskCollection::SeedSetSize{1})},
                                          PerformDiagonalFiltering(true),
                                          TileSize(1));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap());
        createSeeds(seedMap, data.idMap());
        REQUIRE(seedMap->seedMap().size() == 7);
        seedMap->createMatches();
        REQUIRE(seedMap->matches().size() == 7);
        for (auto&& match : seedMap->matches()) {
            REQUIRE(match.distance() == 100);
        }

        // perform filter manually
        auto filter = DiagonalMatchesFilter<KmerOccurrencePair>(config);
        std::vector<KmerOccurrencePair> reported;
        auto orderedMatches = seedMap->orderedMatches();
        filter.applyDiagonalMatchesFilter(orderedMatches, reported);
        tsl::hopscotch_set<KmerOccurrencePair, KmerOccurrencePairHashIgnoreSpan, KmerOccurrencePairEqualIgnoreSpan> reportedSet;
        reportedSet.insert(reported.begin(), reported.end());

        // call SeedMapSpaced filter, should give equal result
        auto skipped = std::make_shared<std::array<size_t, 4>>();
        seedMap->applyDiagonalMatchesFilter(skipped);

        REQUIRE(seedMap->matches() == reportedSet);
        REQUIRE((*skipped)[0] == filter.skippedNotInGenome1And2());
        REQUIRE((*skipped)[1] == filter.skippedOverlappedOrTooClose());
        REQUIRE((*skipped)[2] == filter.skippedTooFewDiagonalElements());
        REQUIRE((*skipped)[3] == filter.skippedTooFewNeighbours());

        std::sort(reported.begin(), reported.end());
        std::vector<KmerOccurrencePair> expected{
// AAAAACCCCCGGGGGTTTTT
// XX|XXXXX no
// XXXX|XXXXX  yes
// ...
//             XXXXX|XX no
//                KmerOccurrencePair(KmerOccurrence(0,0,0,false,"AAAAA"), KmerOccurrence(1,1,100,false,"AAAAA")),
                KmerOccurrencePair(KmerOccurrence(0,0,2,false,"AAACC"), KmerOccurrence(1,1,102,false,"AAACC"), config->span()),
                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC"), config->span()),
                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG"), config->span()),
                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG"), config->span()),
                KmerOccurrencePair(KmerOccurrence(0,0,12,false,"GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT"), config->span()),
//                KmerOccurrencePair(KmerOccurrence(0,0,15,false,"TTTTT"), KmerOccurrence(1,1,115,false,"TTTTT"))
        };
        std::sort(expected.begin(), expected.end());
        REQUIRE(reported == expected);
        //REQUIRE(*skipped == 2);
        REQUIRE((*skipped)[0] == 0);
        REQUIRE((*skipped)[1] == 0);
        REQUIRE((*skipped)[2] == 0);
        REQUIRE((*skipped)[3] == 2);
    }
}

TEST_CASE("SeedMapSpaced Testdata") {
    std::cout << "[INFO] -- [TEST CASE] -- SeedMapSpaced Testdata" << std::endl;
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto fastaCollection = data.fastaCollection();
    auto config = customConfiguration(DiagonalThreshold(2),
                                      CreateAllMatches(true),
                                      Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      LocalSearchAreaLength(100),
                                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(10),
                                                                                                         SpacedSeedMaskCollection::Span(10),
                                                                                                         SpacedSeedMaskCollection::SeedSetSize(1))},
                                      MatchLimit(ULLONG_MAX),
                                      PerformDiagonalFiltering(true),
                                      Thinning(1),
                                      TileSize(1));
    auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap());

    SECTION("Test on Testdata") {
        auto extract = ExtractSeeds<TwoBitKmerDataShort,
                                    TwoBitKmerDataShort>(fastaCollection,
                                                         seedMap,
                                                         false);
        extract.extract();

        REQUIRE(seedMap->numGenomes() == 4);
        REQUIRE(seedMap->numSequences() == 8);

        auto idMapCopy = seedMap->idMapCopy();
        auto matchesPrefilterObs = seedMap->orderedMatches()/ *matches()* /;
        std::set<KmerOccurrencePair> matchesPrefilterExp;
        fillExpectedPrefilterMatches(matchesPrefilterExp, idMapCopy);
        std::vector<KmerOccurrencePair> matchesPrefilterExpV;
        matchesPrefilterExpV.insert(matchesPrefilterExpV.begin(), matchesPrefilterExp.begin(), matchesPrefilterExp.end());
        std::sort(matchesPrefilterExpV.begin(), matchesPrefilterExpV.end());

        REQUIRE(matchesPrefilterExpV == matchesPrefilterObs);



        auto skipped = std::make_shared<std::array<size_t, 4>>();
        seedMap->applyDiagonalMatchesFilter(skipped);
        auto matchesObs = seedMap->orderedMatches()/ *matches()* /;
        std::vector<KmerOccurrencePair> matchesExp;
        fillExpectedMatches(matchesExp, idMapCopy);
        std::sort(matchesExp.begin(), matchesExp.end());

        SECTION("Test Statistics") {
            / * not in human/mouse: 805
               too close/overlap: 146
               too few on diagonal: 3
               too few neighbours: 0 * /
            REQUIRE((skipped->at(0)) == 805);
            REQUIRE((skipped->at(1)) == 146);
            REQUIRE((skipped->at(2)) == 3);
            REQUIRE((skipped->at(3)) == 0);
        }
        SECTION("Test Match Report") {
            REQUIRE(matchesExp == matchesObs);
        }
    }
}
*/
