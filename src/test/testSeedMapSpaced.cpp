#include <bitset>
#include <cstdlib>
#include <climits>
#include <thread>
#include <vector>

#include "catch2/catch.hpp"
#include "prettyprint/prettyprint.hpp"
#include "../ExtractSeeds.h"
#include "../SeedMapSpaced.h"
#include "ConfigurationGenerator.h"
#include "fillSeedMapSpacedTestdata.h"
#include "LoadTestdata.h"

TEST_CASE("SeedMapSpaced") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto config = customConfiguration(Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      SeedSetSize(3),
                                      Span(64),
                                      TileSize(1),
                                      Weight(5));
    auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                  SpacedSeedMaskCollection::Span(config->span()),
                                                                  SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
    auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataLong, TwoBitKmerDataShort>>(config, data.idMap(), masks);

    data.idMap()->querySequenceID("first_sequence", "hg38_orthologs");
    data.idMap()->querySequenceID("second_sequence", "mm10_orthologs");

    SECTION("Test SeedMapSpaced Methods") {
        seedMap->createSeed(TwoBitKmer<TwoBitKmerDataLong>("ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"), KmerOccurrence(0,0,0,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"));
        REQUIRE(seedMap->seedMap().size() == 3); // three different masks, MAY induce the same seed so this may fail from time to time

        seedMap->createSeed(TwoBitKmer<TwoBitKmerDataLong>("ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"), KmerOccurrence(1,1,100,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"));
        seedMap->createMatches();
        REQUIRE(seedMap->matches().size() == 1);    // is a set, match only counted once
        for (auto&& match : seedMap->matches()) {
            REQUIRE(match == KmerOccurrencePair{KmerOccurrence(0,0,0,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"),
                                                KmerOccurrence(1,1,100,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA")});
        }
    }
}



TEST_CASE("SeedMapSpaced Test Data") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto config = customConfiguration(CreateAllMatches(true),
                                      Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      SeedSetSize(2),
                                      Span(20),
                                      TileSize(1),
                                      Weight(10));
    auto masks = std::make_shared<SpacedSeedMaskCollection>(std::vector<std::string>{"10011010011010010011", "01101101001101101000"});
    auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort,
                                                  TwoBitKmerDataShort>>(config, data.idMap(), masks);

    SECTION("Test on Testdata") {
        auto p = (std::thread::hardware_concurrency() > 0)
                ? std::thread::hardware_concurrency()
                : 1;
        if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }
        auto fastaCollection = std::make_shared<FastaCollection>();
        for (auto&& file : inputFiles) { fastaCollection->emplace(FastaRepresentation::genomeFromFilename(file),
                                                                  file); }
        auto extract = ExtractSeeds<TwoBitKmerDataShort,
                                    TwoBitKmerDataShort>(fastaCollection,
                                                         seedMap,
                                                         p, false);
        extract.extract();

        REQUIRE(seedMap->numGenomes() == 4);
        REQUIRE(seedMap->numSequences() == 8);

        auto idMap = seedMap->idMapCopy();
        SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>::SeedMapType seedMapExp;
        fillExpectedSeedMap(seedMapExp, idMap);

        auto & seedMapObs = seedMap->seedMap();
        //SECTION("Expecting Seeds") {
            for (auto&& kmer : seedMapExp) {
                REQUIRE(seedMapObs.find(kmer.first) != seedMapObs.end());
            }
        //}
        //SECTION("Observed Seeds") {
            for (auto&& kmer : seedMapObs) {
                REQUIRE(seedMapExp.find(kmer.first) != seedMapExp.end());
            }
        //}
        //SECTION("Expecting Occs") {
            for (auto&& elem : seedMapExp) {
                for (size_t i = 0; i < masks->size(); ++i) {
                    auto& occs = elem.second.at(i);
                    for (auto&& occ : occs) {
                        REQUIRE(std::find(seedMapObs.at(elem.first).at(i).begin(),
                                          seedMapObs.at(elem.first).at(i).end(),
                                          occ) != seedMapObs.at(elem.first).at(i).end());
                    }
                }
            }
        //}
        //SECTION("Observed Occs") {
            for (auto&& elem : seedMapObs) {
                for (size_t i = 0; i < masks->size(); ++i) {
                    auto& occs = elem.second.at(i);
                    for (auto&& occ : occs) {
                        REQUIRE(std::find(seedMapExp.at(elem.first).at(i).begin(),
                                          seedMapExp.at(elem.first).at(i).end(),
                                          occ) != seedMapExp.at(elem.first).at(i).end());
                    }
                }
            }
        //}

        SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>::MatchesType matchesExp;
        fillExpectedSeedMatches(matchesExp, idMap);
        auto& matchesObs = seedMap->matches();
        //SECTION("Expecting Matches") {
            for (auto&& elem : matchesExp) {
                REQUIRE(matchesObs.find(elem) != matchesObs.end());
            }
        //}
        //SECTION("Observing Matches") {
            for (auto&& elem : matchesObs) {
                REQUIRE(matchesExp.find(elem) != matchesExp.end());
            }
        //}
    }
}
