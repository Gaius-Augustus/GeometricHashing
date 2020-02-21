#include <bitset>
#include <cstdlib>
#include <climits>
#include <thread>
#include <vector>

#include "catch2/catch.hpp"
#include "../ExtractSeeds.h"
#include "../SeedMapContiguous.h"
#include "fillSeedMapContiguousTestdata.h"
#include "LoadTestdata.h"

TEST_CASE("ExactMatches") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto seedMap = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>>(
                       5, genome0, genome1, data.idMap(), ULLONG_MAX, false, false, std::thread::hardware_concurrency()
                   );

    data.idMap()->querySequenceID("first_sequence", "hg38_orthologs");
    data.idMap()->querySequenceID("second_sequence", "mm10_orthologs");

    SECTION("Test ExactMatches Methods") {
        seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(0,0,0,false,"ACGTA"));
        REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).size() == 1);
        REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).size() == 1);
        REQUIRE(seedMap->seedMap().at(TwoBitKmer<TwoBitKmerDataShort>("ACGTA")).at(0).at(0) == KmerOccurrence(0,0,0,false,"ACGTA"));

        seedMap->createSeed(TwoBitKmer<TwoBitKmerDataShort>("ACGTA"), KmerOccurrence(1,1,100,false,"ACGTA"));
        seedMap->createMatches();
        REQUIRE(seedMap->matches().size() == 1);
        REQUIRE(seedMap->matches().at(0) == KmerOccurrencePair{KmerOccurrence(0,0,0,false,"ACGTA"),
                                                               KmerOccurrence(1,1,100,false,"ACGTA")});
    }
}



TEST_CASE("ExactMatches Testdata") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();

    auto kmerMap = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort,
                                                      TwoBitKmerDataShort>>(20, genome0, genome1,
                                                                            data.idMap(),
                                                                            ULLONG_MAX, false,
                                                                            true, std::thread::hardware_concurrency());

    SECTION("Test on Testdata") {
        auto p = (std::thread::hardware_concurrency() > 0)
                ? std::thread::hardware_concurrency()
                : 1;
        if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }
        auto fastaCollection = std::make_shared<FastaCollection>();
        for (auto&& file : inputFiles) { fastaCollection->emplace(FastaRepresentation::genomeFromFilename(file),
                                                                  file); }
        ExtractSeeds<TwoBitKmerDataShort,
                     TwoBitKmerDataShort>(fastaCollection,
                                          kmerMap,
                                          p, false);

        REQUIRE(kmerMap->numGenomes() == 4);
        REQUIRE(kmerMap->numSequences() == 8);

        auto idMap = kmerMap->idMapCopy();
        SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>::SeedMapType kmerMapExp;
        fillExpectedKmerMap(kmerMapExp, idMap);

        auto & kmerMapObs = kmerMap->seedMap();
        for (auto&& kmer : kmerMapExp) {
            REQUIRE(kmerMapObs.find(kmer.first) != kmerMapObs.end());
//            std::cout << "Expected kmer: " << kmer.first << std::endl;
        }
        for (auto&& kmer : kmerMapObs) {
//            std::cout << "Observed kmer: " << kmer.first << std::endl;
            REQUIRE(kmerMapExp.find(kmer.first) != kmerMapExp.end());
        }
        for (auto&& elem : kmerMapExp) {
            auto& occs = elem.second.at(0);
            for (auto&& occ : occs) {
                REQUIRE(std::find(kmerMapObs.at(elem.first).at(0).begin(),
                                  kmerMapObs.at(elem.first).at(0).end(),
                                  occ) != kmerMapObs.at(elem.first).at(0).end());
            }
        }
        for (auto&& elem : kmerMapObs) {
            auto& occs = elem.second.at(0);
            for (auto&& occ : occs) {
                REQUIRE(std::find(kmerMapExp.at(elem.first).at(0).begin(),
                                  kmerMapExp.at(elem.first).at(0).end(),
                                  occ) != kmerMapExp.at(elem.first).at(0).end());
            }
        }

        //kmerMap->createMatches();
        SeedMapContiguous<TwoBitKmerDataShort, TwoBitKmerDataShort>::MatchesType matchesExp;
        fillExpectedMatches(matchesExp, idMap);
        auto& matchesObs = kmerMap->matches();
        REQUIRE(matchesObs.size() == matchesExp.size());
        for (auto&& elem : matchesExp) {
            REQUIRE(std::find(matchesObs.begin(), matchesObs.end(), elem) != matchesObs.end());
        }
        for (auto&& elem : matchesObs) {
            REQUIRE(std::find(matchesExp.begin(), matchesExp.end(), elem) != matchesExp.end());
        }
    }
}
