#include <algorithm>
#include <bitset>
#include <cstdlib>
#include <vector>

#include "catch2/catch.hpp"
#include "../SeedMap.h"
#include "ConfigurationGenerator.h"
#include "LoadTestdata.h"

template<typename SeedMapType>
void createSeed(std::string const & kmer, KmerOccurrence occ,
                std::shared_ptr<SeedMapType> seedMap,
                std::shared_ptr<IdentifierMapping const> idMap,
                std::shared_ptr<Configuration const> config) {
    auto wrapper = [&seedMap](TwoBitKmer<typename SeedMapType::TwoBitSeedDataType_> const & seed, KmerOccurrence const & occurrence, size_t maskIndex) {
        seedMap->addSeed(seed, occurrence, maskIndex);
    };

    ExtractSeeds<typename SeedMapType::TwoBitSeedDataType_> extractor{config->maskCollection(), idMap, config};
    extractor.template createSeed<KmerOccurrence>(kmer, occ, wrapper);
}



TEST_CASE("Test SeedMap 1") {
    std::cout << "[INFO] -- [TEST CASE] -- Test SeedMap 1" << std::endl;
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto config = customConfiguration(Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(16),
                                                                                                         SpacedSeedMaskCollection::Span(20),
                                                                                                         SpacedSeedMaskCollection::SeedSetSize(3))},
                                      NThreads(1),
                                      TileSize(1));
    auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, data.idMap());

    for (auto&& file : inputFiles) {
        REQUIRE(fs::exists(file));
    }

    SECTION("Test SeedMap Methods") {
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



TEST_CASE("Test SeedMap 2") {
    std::cout << "[INFO] -- [TEST CASE] -- Test SeedMap 2" << std::endl;
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto config = customConfiguration(Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{5},
                                                                                                         SpacedSeedMaskCollection::Span{5},
                                                                                                         SpacedSeedMaskCollection::SeedSetSize{1})},
                                      TileSize(1));
    auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, data.idMap());

    data.idMap()->queryGenomeID("third_genome");
    data.idMap()->queryGenomeID("fourth_genome");
    data.idMap()->querySequenceID("first_sequence", "hg38_orthologs");

    SECTION("Test SeedMap Methods") {
        REQUIRE(seedMap->span() == 5);
        REQUIRE(seedMap->weight() == 5);
        REQUIRE(seedMap->seedSetSize() == 1);

        //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(0,0,0,false,"ACGTA"));
        createSeed("ACGTA", KmerOccurrence(0,0,0,false,"ACGTA"), seedMap, data.idMap(), config);
        REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).size() == 1);
        REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).size() == 1);
        REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).at(0) == KmerOccurrence(0,0,0,false,"ACGTA"));

        SECTION("Test missing cleanup mechanisms") {
            //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(0,0,100,false,"ACGTA"));
            createSeed("ACGTA", KmerOccurrence(0,0,100,false,"ACGTA"), seedMap, data.idMap(), config);
            REQUIRE(seedMap->seedMap().size() == 1);

            //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(0,0,0,false,"ACGTC"));
            //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(1,0,0,false,"ACGTC"));
            //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(2,0,0,false,"ACGTC"));
            //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(3,0,0,false,"ACGTC"));
            createSeed("ACGTC", KmerOccurrence(0,0,0,false,"ACGTC"), seedMap, data.idMap(), config);
            createSeed("ACGTC", KmerOccurrence(1,0,0,false,"ACGTC"), seedMap, data.idMap(), config);
            createSeed("ACGTC", KmerOccurrence(2,0,0,false,"ACGTC"), seedMap, data.idMap(), config);
            createSeed("ACGTC", KmerOccurrence(3,0,0,false,"ACGTC"), seedMap, data.idMap(), config);
            REQUIRE(seedMap->seedMap().size() == 2);
        }

        SECTION("Test merging") {
            //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(0,0,0,false,"ACGTC"));
            //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(1,0,0,false,"ACGTC"));
            //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(2,0,0,false,"ACGTC"));
            createSeed("ACGTC", KmerOccurrence(0,0,0,false,"ACGTC"), seedMap, data.idMap(), config);
            createSeed("ACGTC", KmerOccurrence(1,0,0,false,"ACGTC"), seedMap, data.idMap(), config);
            createSeed("ACGTC", KmerOccurrence(2,0,0,false,"ACGTC"), seedMap, data.idMap(), config);

            auto seedMap2 = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, data.idMap());
            //seedMap2->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTG"), KmerOccurrence(1,0,0,false,"ACGTG"));
            //seedMap2->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(0,0,0,false,"ACGTA"));
            //seedMap2->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTC"), KmerOccurrence(3,0,0,false,"ACGTC"));
            createSeed("ACGTG", KmerOccurrence(1,0,0,false,"ACGTG"), seedMap2, data.idMap(), config);
            createSeed("ACGTA", KmerOccurrence(0,0,0,false,"ACGTA"), seedMap2, data.idMap(), config);
            createSeed("ACGTC", KmerOccurrence(3,0,0,false,"ACGTC"), seedMap2, data.idMap(), config);

            seedMap->merge(seedMap2);
            REQUIRE(seedMap->seedMap().size() == 3);
            REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).size() == 1);
            REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).size() == 2);
            REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTC")).size() == 1);
            REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTC")).at(0).size() == 4);
            REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTG")).size() == 1);
            REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTG")).at(0).size() == 1);
            REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTG")).at(0).at(0) == KmerOccurrence(1,0,0,false,"ACGTG"));

            // Test clear
            seedMap->clear();
            REQUIRE(seedMap->size() == 0);
        }
    }
}



TEST_CASE("SeedMap with Spaced") {
    std::cout << "[INFO] -- [TEST CASE] -- SeedMap with Spaced" << std::endl;
    LoadTestdata data;
    auto config = customConfiguration(Genome1(data.genome0()),
                                      Genome2(data.genome1()),
                                      InputFiles(data.inputFiles()),
                                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(5),
                                                                                                         SpacedSeedMaskCollection::Span(64),
                                                                                                         SpacedSeedMaskCollection::SeedSetSize(3))},
                                      TileSize(1));
    auto masks = config->maskCollection();
    auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(config, data.idMap());
    data.idMap()->querySequenceID("first_sequence", "hg38_orthologs");
    REQUIRE(seedMap->seedSetSize() == 3);
    REQUIRE(seedMap->masks().size() == 3);
    //seedMap->createSeed(TwoBitKmer<TwoBitKmerDataLong>("ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"), KmerOccurrence(0,0,0,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"));
    createSeed("ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA", KmerOccurrence(0,0,0,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"), seedMap, data.idMap(), config);
    REQUIRE(seedMap->seedMap().size() == 3);  // MAY fail if more than one seed induce the same sequence
    std::map<size_t, bool> maskUsed;
    for (auto&& seed : seedMap->seedMap()) {
        REQUIRE(seed.first.length() == 5);
        REQUIRE(seed.second.size() == 3);
        size_t lenSum = 0;
        for (size_t i = 0; i < masks->size(); ++i) {
            auto& occVec = seed.second.at(i);
            lenSum += occVec.size();
            if (occVec.size() == 1) {
                maskUsed[i] = true;
                REQUIRE(occVec.at(0) == KmerOccurrence(0,0,0,false,"ACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGAACGTAACGTATAACGA"));
            }
        }
        REQUIRE(lenSum == 1);
    }
    for (size_t i = 0; i < masks->size(); ++i) {
        REQUIRE(maskUsed.at(i));
    }
}
