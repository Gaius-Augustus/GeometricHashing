#include <bitset>
#include <cstdlib>
#include <vector>

#include "catch2/catch.hpp"
#include "../SeedMapGeneral.h"
#include "../SeedMapContiguous.h"
#include "../SeedMapSpaced.h"
#include "ConfigurationGenerator.h"
#include "LoadTestdata.h"

TEST_CASE("Former SeedMapBase") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto config = customConfiguration(Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      NThreads(1),
                                      SeedSetSize(3),
                                      Span(20),
                                      TileSize(1),
                                      Weight(16));
    auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                  SpacedSeedMaskCollection::Span(config->span()),
                                                                  SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
    auto seedMap = std::static_pointer_cast<SeedMapGeneral<TwoBitKmerDataShort, TwoBitKmerDataShort>>(
                       std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks)
                   );

    for (auto&& file : inputFiles) {
        REQUIRE(fs::exists(file));
    }

    SECTION("Test Former SeedMapBase Methods") {
        REQUIRE(seedMap->genome0() == genome0);
        REQUIRE(seedMap->genome1() == genome1);

        auto idMap = seedMap->idMap();
        auto idMapCopy = seedMap->idMapCopy();
        REQUIRE(*idMap == idMapCopy);
        idMapCopy.queryGenomeID("SOME_NEW_GENOME");
        REQUIRE(!(*idMap == idMapCopy));
        REQUIRE(idMap->queryGenomeIDConst(genome0) == 0);
        REQUIRE(idMap->queryGenomeIDConst(genome1) == 1);

        REQUIRE(seedMap->span() == 20);
        REQUIRE(seedMap->numGenomes() == 4);
        REQUIRE(seedMap->numSequences() == 8);
        REQUIRE(seedMap->numGenomes() == idMap->numGenomes());
        REQUIRE(seedMap->numSequences() == idMap->numSequences());
    }
}



TEST_CASE("SeedMapGeneral") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto config = customConfiguration(Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      SeedSetSize(1),
                                      Span(5),
                                      TileSize(1),
                                      Weight(5));
    auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                  SpacedSeedMaskCollection::Span(config->span()),
                                                                  SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
    auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks);
    auto seedMapGeneral = std::static_pointer_cast<SeedMapGeneral<TwoBitKmerDataShort, TwoBitKmerDataShort>>(seedMap);

    data.idMap()->queryGenomeID("third_genome");
    data.idMap()->queryGenomeID("fourth_genome");
    data.idMap()->querySequenceID("first_sequence", "hg38_orthologs");

    SECTION("Test SeedMapGeneral Methods") {
        REQUIRE(seedMapGeneral->span() == 5);
        REQUIRE(seedMapGeneral->weight() == 5);
        REQUIRE(seedMapGeneral->seedSetSize() == 1);

        //seedMapGeneral->updateInvalidSeedStatistics(genome0);
        //REQUIRE(std::get<1>(seedMapGeneral->perGenomeStatistics().at(genome0)) == 1);
        //seedMapGeneral->updateSequenceLengthStatistics(genome0, 42);
        //REQUIRE(std::get<0>(seedMapGeneral->perGenomeStatistics().at(genome0)) == 42);

        seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(0,0,0,false,"ACGTA"));
        REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).size() == 1);
        REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).size() == 1);
        REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).at(0) == KmerOccurrence(0,0,0,false,"ACGTA"));

        SECTION("Test missing cleanup mechanisms") {
            seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(0,0,100,false,"ACGTA"));
            REQUIRE(seedMapGeneral->seedMap().size() == 1);

            seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(0,0,0,false,"ACGTC"));
            seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(1,0,0,false,"ACGTC"));
            seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(2,0,0,false,"ACGTC"));
            seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(3,0,0,false,"ACGTC"));
            REQUIRE(seedMapGeneral->seedMap().size() == 2);
        }

        SECTION("Test merging") {
            seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(0,0,0,false,"ACGTC"));
            seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(1,0,0,false,"ACGTC"));
            seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(2,0,0,false,"ACGTC"));

            auto seedMap2 = seedMap->localCopy();
            seedMap2->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTG"), KmerOccurrence(1,0,0,false,"ACGTG"));
            // after merging, both seeds will be excluded vvv
            seedMap2->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(0,0,0,false,"ACGTA"));
            seedMap2->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(3,0,0,false,"ACGTC"));

            seedMap->merge(seedMap2);
            REQUIRE(seedMapGeneral->seedMap().size() == 3);
            REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).size() == 1);
            REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).size() == 2);
            REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTC")).size() == 1);
            REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTC")).at(0).size() == 4);
            REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTG")).size() == 1);
            REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTG")).at(0).size() == 1);
            REQUIRE(seedMapGeneral->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTG")).at(0).at(0) == KmerOccurrence(1,0,0,false,"ACGTG"));
        }
    }
}

TEST_CASE("SeedMapGeneral with Spaced") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto config = customConfiguration(Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      SeedSetSize(1),
                                      Span(5),
                                      TileSize(1),
                                      Weight(5));
    auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                  SpacedSeedMaskCollection::Span(config->span()),
                                                                  SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
    auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks);
    auto seedMapSpaced = std::static_pointer_cast<SeedMapGeneral<TwoBitKmerDataShort, TwoBitKmerDataShort>>(seedMap);

    data.idMap()->querySequenceID("first_sequence", "hg38_orthologs");

    SECTION("Test SeedMapSpaced Methods") {
        REQUIRE(seedMapSpaced->span() == 5);
        REQUIRE(seedMapSpaced->weight() == 5);
        REQUIRE(seedMapSpaced->seedSetSize() == 1);
        seedMapSpaced->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(0,0,0,false,"ACGTA"));
        REQUIRE(seedMapSpaced->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).size() == 1);
        REQUIRE(seedMapSpaced->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).size() == 1);
        REQUIRE(seedMapSpaced->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).at(0) == KmerOccurrence(0,0,0,false,"ACGTA"));
    }

    SECTION("Test Spaced Seed Methods") {
        auto config2 = customConfiguration(Genome1(data.genome0()),
                                           Genome2(data.genome1()),
                                           InputFiles(data.inputFiles()),
                                           SeedSetSize(3),
                                           Span(64),
                                           TileSize(1),
                                           Weight(5));
        auto masks2 = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config2->weight()),
                                                                       SpacedSeedMaskCollection::Span(config2->span()),
                                                                       SpacedSeedMaskCollection::SeedSetSize(config2->seedSetSize()));
        auto seedMap2 = std::make_shared<SeedMapSpaced<TwoBitKmerDataLong, TwoBitKmerDataShort>>(config2, data.idMap(), masks2);
        auto seedMapSpaced2 = std::static_pointer_cast<SeedMapGeneral<TwoBitKmerDataLong, TwoBitKmerDataShort>>(seedMap2);
        data.idMap()->querySequenceID("first_sequence", "hg38_orthologs");
        REQUIRE(seedMapSpaced2->seedSetSize() == 3);
        REQUIRE(seedMapSpaced2->masks().size() == 3);
        seedMapSpaced2->createSeed(TwoBitKmer<TwoBitKmerDataLong>("ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"), KmerOccurrence(0,0,0,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"));
        REQUIRE(seedMapSpaced2->seedMap().size() == 3);  // MAY fail if more than one seed induce the same sequence
        std::map<size_t, bool> maskUsed;
        for (auto&& seed : seedMapSpaced2->seedMap()) {
            REQUIRE(seed.first.length() == 5);
            REQUIRE(seed.second.size() == 3);
            size_t lenSum = 0;
            for (size_t i = 0; i < masks2->size(); ++i) {
                auto& occVec = seed.second.at(i);
                lenSum += occVec.size();
                if (occVec.size() == 1) {
                    maskUsed[i] = true;
                    REQUIRE(occVec.at(0) == KmerOccurrence(0,0,0,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"));
                }
            }
            REQUIRE(lenSum == 1);
        }
        for (size_t i = 0; i < masks2->size(); ++i) {
            REQUIRE(maskUsed.at(i));
        }
    }
}


TEST_CASE("Match Filter") {
    auto idMap = std::make_shared<IdentifierMapping>("genome0");
    idMap->queryGenomeID("genome1");
    idMap->querySequenceID("seq0", "genome0");
    idMap->querySequenceID("seq1", "genome1");
    SECTION("Two Genomes") {
        // do not create all
        auto seedMap = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>>(4, "genome0", "genome1", idMap,
                                                                                                     2, false, false, std::thread::hardware_concurrency());
        // create all, should not make a difference in this setting
        auto seedMap2 = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>>(4, "genome0", "genome1", idMap,
                                                                                                      2, false, true, std::thread::hardware_concurrency());
        // discard filtered
        auto seedMap3 = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>>(4, "genome0", "genome1", idMap,
                                                                                                      2, true, false, std::thread::hardware_concurrency());
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

        seedMap->createMatches();
        seedMap2->createMatches();
        seedMap3->createMatches();
        REQUIRE(seedMap->matches().size() == 4);
        REQUIRE(seedMap2->matches().size() == 4);
        REQUIRE(seedMap3->matches().size() == 2);

        auto allPossible = std::set<KmerOccurrencePair>{KmerOccurrencePair(o11, o21),
                                                        KmerOccurrencePair(o11, o22),
                                                        KmerOccurrencePair(o12, o21),
                                                        KmerOccurrencePair(o12, o22),
                                                        KmerOccurrencePair(o13, o21),
                                                        KmerOccurrencePair(o13, o22),
                                                        KmerOccurrencePair(o14, o23),
                                                        KmerOccurrencePair(o14, o24)};
        REQUIRE(std::find(seedMap->matches().begin(), seedMap->matches().end(), KmerOccurrencePair(o14, o23)) != seedMap->matches().end()); // find two unfiltered matches
        REQUIRE(std::find(seedMap->matches().begin(), seedMap->matches().end(), KmerOccurrencePair(o14, o24)) != seedMap->matches().end());
        for (auto&& match : seedMap->matches()) {
            REQUIRE(allPossible.find(match) != allPossible.end());  // remaining matches must be from the set of possible matches
        }
        REQUIRE(std::find(seedMap2->matches().begin(), seedMap2->matches().end(), KmerOccurrencePair(o14, o23)) != seedMap2->matches().end());
        REQUIRE(std::find(seedMap2->matches().begin(), seedMap2->matches().end(), KmerOccurrencePair(o14, o24)) != seedMap2->matches().end());
        for (auto&& match : seedMap2->matches()) {
            REQUIRE(allPossible.find(match) != allPossible.end());
        }
        REQUIRE(std::find(seedMap3->matches().begin(), seedMap3->matches().end(), KmerOccurrencePair(o14, o23)) != seedMap3->matches().end());
        REQUIRE(std::find(seedMap3->matches().begin(), seedMap3->matches().end(), KmerOccurrencePair(o14, o24)) != seedMap3->matches().end());
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

        SECTION("Create All") {
            // sample
            auto seedMap = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>>(4, "genome0", "genome1", idMap,
                                                                                                         3, false, true, std::thread::hardware_concurrency());
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // sample 3 from here
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o31, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o32, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o41, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o42, 0);

            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o13, 0); // make all
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o43, 0);
            seedMap->createMatches();

            // discard
            auto seedMap2 = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>>(4, "genome0", "genome1", idMap,
                                                                                                          3, true, true, std::thread::hardware_concurrency());
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // discard
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o31, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o32, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o41, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o42, 0);

            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o13, 0); // make all
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o43, 0);
            seedMap2->createMatches();

            auto allPossible = std::set<KmerOccurrencePair>{KmerOccurrencePair(o11, o21),
                                                            KmerOccurrencePair(o11, o22),
                                                            KmerOccurrencePair(o12, o21),
                                                            KmerOccurrencePair(o12, o22),
                                                            KmerOccurrencePair(o11, o31),
                                                            KmerOccurrencePair(o11, o32),
                                                            KmerOccurrencePair(o12, o31),
                                                            KmerOccurrencePair(o12, o32),
                                                            KmerOccurrencePair(o11, o41),
                                                            KmerOccurrencePair(o11, o42),
                                                            KmerOccurrencePair(o12, o41),
                                                            KmerOccurrencePair(o12, o42),

                                                            KmerOccurrencePair(o21, o31),
                                                            KmerOccurrencePair(o21, o32),
                                                            KmerOccurrencePair(o22, o31),
                                                            KmerOccurrencePair(o22, o32),
                                                            KmerOccurrencePair(o21, o41),
                                                            KmerOccurrencePair(o21, o42),
                                                            KmerOccurrencePair(o22, o41),
                                                            KmerOccurrencePair(o22, o42),


                                                            KmerOccurrencePair(o31, o41),
                                                            KmerOccurrencePair(o31, o42),
                                                            KmerOccurrencePair(o32, o41),
                                                            KmerOccurrencePair(o32, o42),


                                                            KmerOccurrencePair(o13, o23),
                                                            KmerOccurrencePair(o13, o43),
                                                            KmerOccurrencePair(o23, o43)};

            REQUIRE(seedMap->matches().size() == 6);
            REQUIRE(std::find(seedMap->matches().begin(), seedMap->matches().end(), KmerOccurrencePair(o13, o23)) != seedMap->matches().end()); // find unfiltered matches
            REQUIRE(std::find(seedMap->matches().begin(), seedMap->matches().end(), KmerOccurrencePair(o13, o43)) != seedMap->matches().end());
            REQUIRE(std::find(seedMap->matches().begin(), seedMap->matches().end(), KmerOccurrencePair(o23, o43)) != seedMap->matches().end());
            for (auto&& match : seedMap->matches()) {
                REQUIRE(allPossible.find(match) != allPossible.end());  // remaining matches must be from the set of possible matches
            }
            REQUIRE(seedMap2->matches().size() == 3);
            REQUIRE(std::find(seedMap2->matches().begin(), seedMap2->matches().end(), KmerOccurrencePair(o13, o23)) != seedMap2->matches().end()); // find unfiltered matches
            REQUIRE(std::find(seedMap2->matches().begin(), seedMap2->matches().end(), KmerOccurrencePair(o13, o43)) != seedMap2->matches().end());
            REQUIRE(std::find(seedMap2->matches().begin(), seedMap2->matches().end(), KmerOccurrencePair(o23, o43)) != seedMap2->matches().end());

        }
        SECTION("Not Create All") {
            auto seedMap = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>>(4, "genome0", "genome1", idMap,
                                                                                                         3, false, false, std::thread::hardware_concurrency());
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // only four, sample 3 from here
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o31, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o32, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o41, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o42, 0);

            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o13, 0); // make only o13--o23
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
            seedMap->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o43, 0);
            seedMap->createMatches();

            // discard
            auto seedMap2 = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>>(4, "genome0", "genome1", idMap,
                                                                                                          3, true, false, std::thread::hardware_concurrency());
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o11, 0); // discard
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o12, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o21, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o22, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o31, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o32, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o41, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGT"), o42, 0);

            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o13, 0); // make only o13--o23
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o23, 0);
            seedMap2->addSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGA"), o43, 0);
            seedMap2->createMatches();

            auto allPossible = std::set<KmerOccurrencePair>{KmerOccurrencePair(o11, o21),
                                                            KmerOccurrencePair(o11, o22),
                                                            KmerOccurrencePair(o12, o21),
                                                            KmerOccurrencePair(o12, o22),

                                                            KmerOccurrencePair(o13, o23)};

            REQUIRE(seedMap->matches().size() == 4);
            REQUIRE(std::find(seedMap->matches().begin(), seedMap->matches().end(), KmerOccurrencePair(o13, o23)) != seedMap->matches().end()); // find unfiltered match
            for (auto&& match : seedMap->matches()) {
                REQUIRE(allPossible.find(match) != allPossible.end());  // remaining matches must be from the set of possible matches
            }
            REQUIRE(seedMap2->matches().size() == 1);
            REQUIRE(std::find(seedMap2->matches().begin(), seedMap2->matches().end(), KmerOccurrencePair(o13, o23)) != seedMap2->matches().end()); // find unfiltered match
        }
    }
}
