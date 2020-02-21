#include <bitset>
#include <cstdlib>
#include <vector>

#include "catch2/catch.hpp"
#include "../SeedMapBase.h"
#include "../SeedMapSpaced.h"
#include "LoadTestdata.h"

TEST_CASE("SeedMapBase") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                  //ArtificialSequenceLengths(Configuration::ArtificialSequenceLengthsMapType()),
                                                  ArtificialSequenceSizeFactor(1),
                                                  CreateAllMatches(false),
                                                  CubeScoreMu(3),
                                                  CubeScoreThreshold(ULLONG_MAX),
                                                  DiagonalThreshold(0),
                                                  DynamicArtificialSequences(false),
                                                  Genome1(data.genome0()),
                                                  Genome2(data.genome1()),
                                                  GlobalOccurrenceThreshold(ULLONG_MAX),
                                                  InputFiles(data.inputFiles()),
                                                  LinkLimit(ULLONG_MAX),
                                                  LocalSearchAreaLength(0),
                                                  LocalOccurrenceThreshold(1),
                                                  Masks(Configuration::MasksType()),
                                                  MinMatchDistance(0),
                                                  NoProgressbar(false),
                                                  NThreads(1),
                                                  OccurrencePerGenomeMax(0),
                                                  OccurrencePerGenomeMin(0),
                                                  OptimalSeed(false),
                                                  OrthologTuples(Configuration::OrthologTupleListType()),
                                                  OrthologTupleListFile(""),
                                                  Output(""),
                                                  OutputArtificialSequences(""),
                                                  OutputRunInformation(""),
                                                  PerformDiagonalFiltering(false),
                                                  PerformGeometricHashing(false),
                                                  SeedSetSize(3),
                                                  Span(20),
                                                  TileSize(1),
                                                  Weight(16));
    auto masks = std::make_shared<SpacedSeedMaskCollection const>(config->k(), config->l(), config->m());
    auto seedMap = std::static_pointer_cast<SeedMapBase<TwoBitKmerDataShort, TwoBitKmerDataShort>>(
                       std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks)
                   );

    for (auto&& file : inputFiles) {
        REQUIRE(fs::exists(file));
    }

    SECTION("Test SeedMapBase Methods") {
        REQUIRE(seedMap->genome0() == genome0);
        REQUIRE(seedMap->genome1() == genome1);

        auto idMap = seedMap->idMap();
        auto idMapCopy = seedMap->idMapCopy();
        REQUIRE(*idMap == idMapCopy);
        idMapCopy.queryGenomeID("SOME_NEW_GENOME");
        REQUIRE(!(*idMap == idMapCopy));
        REQUIRE(idMap->queryGenomeID(genome0) == 0);
        REQUIRE(idMap->queryGenomeID(genome1) == 1);

        REQUIRE(seedMap->k() == 20);
        REQUIRE(seedMap->numGenomes() == 2);
        REQUIRE(seedMap->numSequences() == 0);
    }
}
