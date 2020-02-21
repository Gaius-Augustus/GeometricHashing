#include <bitset>
#include <cstdlib>
#include <climits>
#include <thread>
#include <vector>

#include "catch2/catch.hpp"
#include "prettyprint/prettyprint.hpp"
#include "../ExtractSeeds.h"
#include "../SeedMapSpaced.h"
#include "fillDiagonalMatchesFilterTestdata.h"
#include "LoadTestdata.h"

void createSeeds(std::shared_ptr<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>> seedMap) {
    seedMap->idMap()->querySequenceID("first_sequence", "hg38_orthologs");
    seedMap->idMap()->querySequenceID("second_sequence", "mm10_orthologs");

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
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();

    SECTION("Test Chunk Creation") {
        auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                      ArtificialSequenceSizeFactor(1),
                                                      CreateAllMatches(false),
                                                      CubeScoreMu(3),
                                                      CubeScoreThreshold(ULLONG_MAX),
                                                      DiagonalThreshold(1),
                                                      DynamicArtificialSequences(false),
                                                      Genome1("genome0"),
                                                      Genome2("genome1"),
                                                      InputFiles(data.inputFiles()),
                                                      LocalSearchAreaLength(1),
                                                      Masks(Configuration::MasksType()),
                                                      MatchLimit(ULLONG_MAX),
                                                      MatchLimitDiscardSeeds(false),
                                                      MinMatchDistance(1),
                                                      NoProgressbar(false),
                                                      NThreads(4),  // needs to be fix in order for chunk algorithm to work predictably
                                                      OccurrencePerGenomeMax(1),
                                                      OccurrencePerGenomeMin(1),
                                                      OptimalSeed(false),
                                                      Output(""),
                                                      OutputArtificialSequences(""),
                                                      OutputRunInformation(""),
                                                      PerformDiagonalFiltering(true),
                                                      PerformGeometricHashing(false),
                                                      SeedSetSize(1),
                                                      Span(5),
                                                      TileSize(1),
                                                      Weight(5));
        auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                      SpacedSeedMaskCollection::Span(config->span()),
                                                                      SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
        auto idMap = std::make_shared<IdentifierMapping>("genome0");
        idMap->queryGenomeID("genome1");
        idMap->querySequenceID("seq0", "genome0");
        idMap->querySequenceID("seq1", "genome1");
        auto seedMap = SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>(config, idMap, masks);
        auto filter = DiagonalMatchesFilter<KmerOccurrencePair>(config);
        SECTION("Test One Match Chunk Creation") {
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> chunk;
            chunk.emplace(KmerOccurrence(0,0,0,false,"ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"));
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> matchSet;
            matchSet.insert(chunk.begin(), chunk.end());
            seedMap.appendMatches(matchSet);
            auto orderedMatches = seedMap.orderedMatches();
            auto chunks = filter.matchesChunks(orderedMatches);
            REQUIRE(chunks.size() == 1);
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reconstructChunk;
            reconstructChunk.insert(chunks.at(0).first, chunks.at(0).second);
            REQUIRE(reconstructChunk == chunk);
        }
        SECTION("Test Two Match Chunk Creation") {
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> chunk;
            chunk.emplace(KmerOccurrence(0,0,0,false,"ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"));
            chunk.emplace(KmerOccurrence(0,0,1,false,"ACGTA"), KmerOccurrence(1,1,101,false,"ACGTA"));
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> matchSet;
            matchSet.insert(chunk.begin(), chunk.end());
            seedMap.appendMatches(matchSet);
            auto orderedMatches = seedMap.orderedMatches();
            auto chunks = filter.matchesChunks(orderedMatches);
            REQUIRE(chunks.size() == 1);
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reconstructChunk;
            reconstructChunk.insert(chunks.at(0).first, chunks.at(0).second);
            REQUIRE(reconstructChunk == chunk);
        }
        SECTION("Test Four Match Chunk Creation (= nThreads_)") {
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> chunk;
            chunk.emplace(KmerOccurrence(0,0,0,false,"ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"));    // distance: +100
            chunk.emplace(KmerOccurrence(0,0,1,false,"ACGTA"), KmerOccurrence(1,1,101,false,"ACGTA"));    // distance: +100
            chunk.emplace(KmerOccurrence(0,0,2,false,"ACGTA"), KmerOccurrence(1,1,102,false,"ACGTA"));    // distance: +100
            chunk.emplace(KmerOccurrence(0,0,3,false,"ACGTA"), KmerOccurrence(1,1,103,false,"ACGTA"));    // distance: +100, all same diag, so one chunk desired
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> matchSet;
            matchSet.insert(chunk.begin(), chunk.end());
            seedMap.appendMatches(matchSet);
            auto orderedMatches = seedMap.orderedMatches();
            auto chunks = filter.matchesChunks(orderedMatches);
            REQUIRE(chunks.size() == 1);
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reconstructChunk;
            reconstructChunk.insert(chunks.at(0).first, chunks.at(0).second);
            REQUIRE(reconstructChunk == chunk);
        }
        SECTION("Test Correct Chunk Creation") {
            // nThreads = 4, wholeSet size = 16, so equal sized would be 4 pairs per chunk
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> chunk0;
            chunk0.emplace(KmerOccurrence(0,0,0,false,"ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"));    // distance: +100
            chunk0.emplace(KmerOccurrence(0,0,1,false,"ACGTA"), KmerOccurrence(1,1,101,false,"ACGTA"));    // distance: +100
            chunk0.emplace(KmerOccurrence(0,0,2,false,"ACGTA"), KmerOccurrence(1,1,102,false,"ACGTA"));    // distance: +100
            chunk0.emplace(KmerOccurrence(0,0,10,false,"ACGTA"), KmerOccurrence(1,1,210,false,"ACGTA"));    // distance: +200 // first equal chunk last
            chunk0.emplace(KmerOccurrence(0,0,11,false,"ACGTA"), KmerOccurrence(1,1,211,false,"ACGTA"));    // distance: +200 // desired first chunk last

            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> chunk1;
            chunk1.emplace(KmerOccurrence(0,0,20,false,"ACGTA"), KmerOccurrence(1,1,320,false,"ACGTA"));    // distance: +300
            chunk1.emplace(KmerOccurrence(0,0,21,false,"ACGTA"), KmerOccurrence(1,1,321,false,"ACGTA"));    // distance: +300
            chunk1.emplace(KmerOccurrence(0,0,22,false,"ACGTA"), KmerOccurrence(1,1,322,false,"ACGTA"));    // distance: +300
            chunk1.emplace(KmerOccurrence(0,0,30,false,"ACGTA"), KmerOccurrence(1,1,430,false,"ACGTA"));    // distance: +400
            chunk1.emplace(KmerOccurrence(0,0,31,false,"ACGTA"), KmerOccurrence(1,1,431,false,"ACGTA"));    // distance: +400
            chunk1.emplace(KmerOccurrence(0,0,32,false,"ACGTA"), KmerOccurrence(1,1,432,false,"ACGTA"));    // distance: +400
            chunk1.emplace(KmerOccurrence(0,0,33,false,"ACGTA"), KmerOccurrence(1,1,433,false,"ACGTA"));    // distance: +400
            chunk1.emplace(KmerOccurrence(0,0,34,false,"ACGTA"), KmerOccurrence(1,1,434,false,"ACGTA"));    // distance: +400 // desired second chunk last

            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> chunk2;
            chunk2.emplace(KmerOccurrence(0,0,40,false,"ACGTA"), KmerOccurrence(1,1,540,false,"ACGTA"));    // distance: +500
            chunk2.emplace(KmerOccurrence(0,0,41,false,"ACGTA"), KmerOccurrence(1,1,541,false,"ACGTA"));    // distance: +500
            chunk2.emplace(KmerOccurrence(0,0,42,false,"ACGTA"), KmerOccurrence(1,1,542,false,"ACGTA"));    // distance: +500

            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> wholeSet;
            wholeSet.insert(chunk0.begin(), chunk0.end());
            wholeSet.insert(chunk1.begin(), chunk1.end());
            wholeSet.insert(chunk2.begin(), chunk2.end());

            seedMap.appendMatches(wholeSet);

            auto orderedMatches = seedMap.orderedMatches();
            auto chunks = filter.matchesChunks(orderedMatches);
            REQUIRE(chunks.size() == 3);
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> builtChunk0;
            builtChunk0.insert(chunks.at(0).first, chunks.at(0).second);
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> builtChunk1;
            builtChunk1.insert(chunks.at(1).first, chunks.at(1).second);
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> builtChunk2;
            builtChunk2.insert(chunks.at(2).first, chunks.at(2).second);
            REQUIRE(builtChunk0 == chunk0);
            REQUIRE(builtChunk1 == chunk1);
            REQUIRE(builtChunk2 == chunk2);
            /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reconstructChunk;
            reconstructChunk.insert(chunks.at(0).first, chunks.at(0).second);
            reconstructChunk.insert(chunks.at(2).first, chunks.at(1).second);
            reconstructChunk.insert(chunks.at(1).first, chunks.at(2).second);
            REQUIRE(reconstructChunk == wholeSet);
        }
    }

    SECTION("Test SeedMapSpaced Methods 1") {
        auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                      ArtificialSequenceSizeFactor(1),
                                                      CreateAllMatches(false),
                                                      CubeScoreMu(3),
                                                      CubeScoreThreshold(ULLONG_MAX),
                                                      DiagonalThreshold(3),
                                                      DynamicArtificialSequences(false),
                                                      Genome1(data.genome0()),
                                                      Genome2(data.genome1()),
                                                      InputFiles(data.inputFiles()),
                                                      LocalSearchAreaLength(15),
                                                      Masks(Configuration::MasksType()),
                                                      MatchLimit(ULLONG_MAX),
                                                      MatchLimitDiscardSeeds(false),
                                                      MinMatchDistance(0),
                                                      NoProgressbar(false),
                                                      NThreads(std::thread::hardware_concurrency()),
                                                      OccurrencePerGenomeMax(1),
                                                      OccurrencePerGenomeMin(1),
                                                      OptimalSeed(false),
                                                      Output(""),
                                                      OutputArtificialSequences(""),
                                                      OutputRunInformation(""),
                                                      PerformDiagonalFiltering(true),
                                                      PerformGeometricHashing(false),
                                                      SeedSetSize(1),
                                                      Span(5),
                                                      TileSize(1),
                                                      Weight(5));
        auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                      SpacedSeedMaskCollection::Span(config->span()),
                                                                      SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks);
        createSeeds(seedMap);
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
        /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reportedSet;
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
                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC")),
//                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG")),
                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG")),
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
        auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                      ArtificialSequenceSizeFactor(1),
                                                      CreateAllMatches(false),
                                                      CubeScoreMu(3),
                                                      CubeScoreThreshold(ULLONG_MAX),
                                                      DiagonalThreshold(3),
                                                      DynamicArtificialSequences(false),
                                                      Genome1(data.genome0()),
                                                      Genome2(data.genome1()),
                                                      InputFiles(data.inputFiles()),
                                                      LocalSearchAreaLength(24),
                                                      Masks(Configuration::MasksType()),
                                                      MatchLimit(ULLONG_MAX),
                                                      MatchLimitDiscardSeeds(false),
                                                      MinMatchDistance(0),
                                                      NoProgressbar(false),
                                                      NThreads(std::thread::hardware_concurrency()),
                                                      OccurrencePerGenomeMax(1),
                                                      OccurrencePerGenomeMin(1),
                                                      OptimalSeed(false),
                                                      Output(""),
                                                      OutputArtificialSequences(""),
                                                      OutputRunInformation(""),
                                                      PerformDiagonalFiltering(true),
                                                      PerformGeometricHashing(false),
                                                      SeedSetSize(1),
                                                      Span(5),
                                                      TileSize(1),
                                                      Weight(5));
        auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                      SpacedSeedMaskCollection::Span(config->span()),
                                                                      SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks);
        createSeeds(seedMap);
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
        filter.applyDiagonalMatchesFilter(orderedMatches/*matches()*/, reported);
        /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reportedSet;
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
                KmerOccurrencePair(KmerOccurrence(0,0,0,false,"AAAAA"), KmerOccurrence(1,1,100,false,"AAAAA")),
//                KmerOccurrencePair(KmerOccurrence(0,0,2,false,"AAACC"), KmerOccurrence(1,1,102,false,"AAACC")),
                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC")),
//                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG")),
                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG")),
//                KmerOccurrencePair(KmerOccurrence(0,0,12,false,"GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT")),
                KmerOccurrencePair(KmerOccurrence(0,0,15,false,"TTTTT"), KmerOccurrence(1,1,115,false,"TTTTT"))
        };
        std::sort(expected.begin(), expected.end());
        REQUIRE(reported == expected);
        REQUIRE((*skipped)[0] == 0);
        REQUIRE((*skipped)[1] == 3);
        REQUIRE((*skipped)[2] == 0);
        REQUIRE((*skipped)[3] == 0);
    }



    SECTION("Test SeedMapSpaced Methods 3") {
        auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                      ArtificialSequenceSizeFactor(1),
                                                      CreateAllMatches(false),
                                                      CubeScoreMu(3),
                                                      CubeScoreThreshold(ULLONG_MAX),
                                                      DiagonalThreshold(3),
                                                      DynamicArtificialSequences(false),
                                                      Genome1(data.genome0()),
                                                      Genome2(data.genome1()),
                                                      InputFiles(data.inputFiles()),
                                                      LocalSearchAreaLength(24),
                                                      Masks(Configuration::MasksType()),
                                                      MatchLimit(ULLONG_MAX),
                                                      MatchLimitDiscardSeeds(false),
                                                      MinMatchDistance(1),
                                                      NoProgressbar(false),
                                                      NThreads(std::thread::hardware_concurrency()),
                                                      OccurrencePerGenomeMax(1),
                                                      OccurrencePerGenomeMin(1),
                                                      OptimalSeed(false),
                                                      Output(""),
                                                      OutputArtificialSequences(""),
                                                      OutputRunInformation(""),
                                                      PerformDiagonalFiltering(true),
                                                      PerformGeometricHashing(false),
                                                      SeedSetSize(1),
                                                      Span(5),
                                                      TileSize(1),
                                                      Weight(5));
        auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                      SpacedSeedMaskCollection::Span(config->span()),
                                                                      SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks);
        createSeeds(seedMap);
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
        /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reportedSet;
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
                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG")),
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
        auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                      ArtificialSequenceSizeFactor(1),
                                                      CreateAllMatches(false),
                                                      CubeScoreMu(3),
                                                      CubeScoreThreshold(ULLONG_MAX),
                                                      DiagonalThreshold(2),
                                                      DynamicArtificialSequences(false),
                                                      Genome1(data.genome0()),
                                                      Genome2(data.genome1()),
                                                      InputFiles(data.inputFiles()),
                                                      LocalSearchAreaLength(24),
                                                      Masks(Configuration::MasksType()),
                                                      MatchLimit(ULLONG_MAX),
                                                      MatchLimitDiscardSeeds(false),
                                                      MinMatchDistance(1),
                                                      NoProgressbar(false),
                                                      NThreads(std::thread::hardware_concurrency()),
                                                      OccurrencePerGenomeMax(1),
                                                      OccurrencePerGenomeMin(1),
                                                      OptimalSeed(false),
                                                      Output(""),
                                                      OutputArtificialSequences(""),
                                                      OutputRunInformation(""),
                                                      PerformDiagonalFiltering(true),
                                                      PerformGeometricHashing(false),
                                                      SeedSetSize(1),
                                                      Span(5),
                                                      TileSize(1),
                                                      Weight(5));
        auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                      SpacedSeedMaskCollection::Span(config->span()),
                                                                      SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks);
        createSeeds(seedMap);
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
        /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reportedSet;
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
                KmerOccurrencePair(KmerOccurrence(0,0,0,false,"AAAAA"), KmerOccurrence(1,1,100,false,"AAAAA")),
//                KmerOccurrencePair(KmerOccurrence(0,0,2,false,"AAACC"), KmerOccurrence(1,1,102,false,"AAACC")),
//                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC")),
                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG")),
//                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG")),
//                KmerOccurrencePair(KmerOccurrence(0,0,12,false,"GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT")),
                KmerOccurrencePair(KmerOccurrence(0,0,15,false,"TTTTT"), KmerOccurrence(1,1,115,false,"TTTTT"))
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
        auto config = std::make_shared<Configuration>(AllowOverlap(true),
                                                      ArtificialSequenceSizeFactor(1),
                                                      CreateAllMatches(false),
                                                      CubeScoreMu(3),
                                                      CubeScoreThreshold(ULLONG_MAX),
                                                      DiagonalThreshold(2),
                                                      DynamicArtificialSequences(false),
                                                      Genome1(data.genome0()),
                                                      Genome2(data.genome1()),
                                                      InputFiles(data.inputFiles()),
                                                      LocalSearchAreaLength(10),
                                                      Masks(Configuration::MasksType()),
                                                      MatchLimit(ULLONG_MAX),
                                                      MatchLimitDiscardSeeds(false),
                                                      MinMatchDistance(0),
                                                      NoProgressbar(false),
                                                      NThreads(std::thread::hardware_concurrency()),
                                                      OccurrencePerGenomeMax(1),
                                                      OccurrencePerGenomeMin(1),
                                                      OptimalSeed(false),
                                                      Output(""),
                                                      OutputArtificialSequences(""),
                                                      OutputRunInformation(""),
                                                      PerformDiagonalFiltering(true),
                                                      PerformGeometricHashing(false),
                                                      SeedSetSize(1),
                                                      Span(5),
                                                      TileSize(1),
                                                      Weight(5));
        auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                      SpacedSeedMaskCollection::Span(config->span()),
                                                                      SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
        auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks);
        createSeeds(seedMap);
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
        /*std::set<KmerOccurrencePair>*/std::unordered_set<KmerOccurrencePair, KmerOccurrencePairHash> reportedSet;
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
                KmerOccurrencePair(KmerOccurrence(0,0,2,false,"AAACC"), KmerOccurrence(1,1,102,false,"AAACC")),
                KmerOccurrencePair(KmerOccurrence(0,0,5,false,"CCCCC"), KmerOccurrence(1,1,105,false,"CCCCC")),
                KmerOccurrencePair(KmerOccurrence(0,0,7,false,"CCCGG"), KmerOccurrence(1,1,107,false,"CCCGG")),
                KmerOccurrencePair(KmerOccurrence(0,0,10,false,"GGGGG"), KmerOccurrence(1,1,110,false,"GGGGG")),
                KmerOccurrencePair(KmerOccurrence(0,0,12,false,"GGGTT"), KmerOccurrence(1,1,112,false,"GGGTT")),
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
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto fastaCollection = data.fastaCollection();
    auto config = std::make_shared<Configuration>(AllowOverlap(false),
                                                  ArtificialSequenceSizeFactor(1),
                                                  CreateAllMatches(true),
                                                  CubeScoreMu(3),
                                                  CubeScoreThreshold(ULLONG_MAX),
                                                  DiagonalThreshold(2.),
                                                  DynamicArtificialSequences(false),
                                                  Genome1(data.genome0()),
                                                  Genome2(data.genome1()),
                                                  InputFiles(data.inputFiles()),
                                                  LocalSearchAreaLength(100),
                                                  Masks(Configuration::MasksType()),
                                                  MatchLimit(ULLONG_MAX),
                                                  MatchLimitDiscardSeeds(false),
                                                  MinMatchDistance(0),
                                                  NoProgressbar(false),
                                                  NThreads(std::thread::hardware_concurrency()),
                                                  OccurrencePerGenomeMax(1),
                                                  OccurrencePerGenomeMin(1),
                                                  OptimalSeed(false),
                                                  Output(""),
                                                  OutputArtificialSequences(""),
                                                  OutputRunInformation(""),
                                                  PerformDiagonalFiltering(true),
                                                  PerformGeometricHashing(false),
                                                  SeedSetSize(1),
                                                  Span(10),
                                                  TileSize(1),
                                                  Weight(10));
    auto masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                  SpacedSeedMaskCollection::Span(config->span()),
                                                                  SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
    auto seedMap = std::make_shared<SeedMapSpaced<TwoBitKmerDataShort, TwoBitKmerDataShort>>(config, data.idMap(), masks);

    SECTION("Test on Testdata") {
        auto p = (std::thread::hardware_concurrency() > 0)
                ? std::thread::hardware_concurrency()
                : 1;
        if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }
        ExtractSeeds<TwoBitKmerDataShort,
                     TwoBitKmerDataShort>(fastaCollection,
                                          seedMap,
                                          p, false);

        REQUIRE(seedMap->numGenomes() == 4);
        REQUIRE(seedMap->numSequences() == 8);

        auto idMapCopy = seedMap->idMapCopy();
        auto matchesPrefilterObs = seedMap->orderedMatches()/*matches()*/;
        std::set<KmerOccurrencePair> matchesPrefilterExp;
        fillExpectedPrefilterMatches(matchesPrefilterExp, idMapCopy);

        REQUIRE(matchesPrefilterExp == matchesPrefilterObs);



        auto skipped = std::make_shared<std::array<size_t, 4>>();
        seedMap->applyDiagonalMatchesFilter(skipped);
        auto matchesObs = seedMap->orderedMatches()/*matches()*/;
        std::vector<KmerOccurrencePair> matchesExpVector;
        fillExpectedMatches(matchesExpVector, idMapCopy);

        std::set<KmerOccurrencePair> matchesExp;
        matchesExp.insert(matchesExpVector.begin(), matchesExpVector.end());
        //std::sort(matchesExp.begin(), matchesExp.end());
        //std::sort(matchesObs.begin(), matchesObs.end());

        SECTION("Test Statistics") {
            /* not in human/mouse: 805
               too close/overlap: 146
               too few on diagonal: 3
               too few neighbours: 0 */
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

