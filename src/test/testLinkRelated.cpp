#include <experimental/filesystem>
#include <thread>
#include <fstream>
#include <memory>

#include "catch2/catch.hpp"
#include "../ExtractSeeds.h"
#include "../IdentifierMapping.h"
#include "../Link.h"
#include "../Linkset.h"
#include "../SeedMapContiguous.h"
#include "fillLinksetTestdata.h"
#include "GeometricHashingTestdata.h"
#include "LoadTestdata.h"
//#include "manualTestset.h"

/* ~~~~~~~~~~~~
 * TESTING Link
 * ~~~~~~~~~~~~ */
TEST_CASE("Link") {
    Link link;

    SECTION("Test Link data") {
        auto thrown = false;
        try {
            link.first();
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);
        thrown = false;
        try {
            link.second();
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);

        link.insertOccurrence(2,2,100,false,"ACGT"); // 3rd after next two insertions

        thrown = false;
        try {
            link.second();
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);

        link.insertOccurrence(1,1,0,true,"ACGT");    // 2nd after next insertion
        link.insertOccurrence(0,0,0,false,"ACGT");   // 1st after insertion

        REQUIRE(link.dimensionality() == 3);
        REQUIRE(link.genome(0) == 0);
        REQUIRE(link.genome(1) == 1);
        REQUIRE(link.genome(2) == 2);
        REQUIRE(!link.reverse(0));
        REQUIRE(link.reverse(1));
        REQUIRE(!link.reverse(2));
        REQUIRE(link.sequence(0) == 0);
        REQUIRE(link.sequence(1) == 1);
        REQUIRE(link.sequence(2) == 2);
        REQUIRE(link.position(0) == 0);
        REQUIRE(link.position(1) == 0);
        REQUIRE(link.position(2) == 100);
        REQUIRE(link.first() == KmerOccurrence{0,0,0,false,"ACGT"});
        REQUIRE(link.second() == KmerOccurrence{1,1,0,true,"ACGT"});

        thrown = false;
        try {
            link.insertOccurrence(1,3,0,true,"ACGT"); // throw: genome 1 twice
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);

        // Construct from vector
        std::vector<KmerOccurrence> vec;
        vec.emplace_back(2,2,100,false,"ACGT"); // 3rd after next two insertions
        vec.emplace_back(1,1,0,true,"ACGT");    // 2nd after next insertion
        vec.emplace_back(0,0,0,false,"ACGT");   // 1st after insertion
        vec.emplace_back(1,1,0,true,"ACGT");    // leads to throw

        thrown = false;
        try {
            auto link2 = Link(vec);
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);

        vec.pop_back();
        auto link2 = Link(vec);
        REQUIRE(link2.dimensionality() == 3);
        REQUIRE(link2.genome(0) == 0);
        REQUIRE(link2.genome(1) == 1);
        REQUIRE(link2.genome(2) == 2);
        REQUIRE(!link2.reverse(0));
        REQUIRE(link2.reverse(1));
        REQUIRE(!link2.reverse(2));
        REQUIRE(link2.sequence(0) == 0);
        REQUIRE(link2.sequence(1) == 1);
        REQUIRE(link2.sequence(2) == 2);

        REQUIRE(link2.sameDistance(link));
    }

    SECTION("Test Link comparison") {
        Link l1;
        Link l2;
        Link l3;
        Link l4;
        Link l5;
        Link l6;
        Link l7;

        l1.insertOccurrence(0, 0, 0, false, "ACGT");
        l1.insertOccurrence(1, 1, 100, false, "ACGT");
        l1.insertOccurrence(2, 1, 100, true, "ACGT");
        l1.insertOccurrence(3, 0, 500, false, "ACGT");

        l2.insertOccurrence(0, 0, 0, false, "ACGT");  // == l1
        l2.insertOccurrence(1, 1, 100, false, "ACGT");
        l2.insertOccurrence(2, 1, 100, true, "ACGT");
        l2.insertOccurrence(3, 0, 500, false, "ACGT");

        l3.insertOccurrence(0, 0, 0, false, "ACGT");  // < l1
        l3.insertOccurrence(1, 1, 100, false, "ACGT");
        l3.insertOccurrence(2, 1, 100, true, "ACGT");

        l4.insertOccurrence(0, 0, 0, false, "ACGT"); // > l1
        l4.insertOccurrence(1, 1, 100, false, "ACGT");
        l4.insertOccurrence(3, 0, 500, true, "ACGT");

        l5.insertOccurrence(0, 0, 0, false, "ACGT");  // == l1
        l5.insertOccurrence(1, 1, 100, false, "ACGT");
        l5.insertOccurrence(2, 1, 100, true, "ACGT");
        l5.insertOccurrence(3, 0, 500, false, "ACGT");

        l6.insertOccurrence(0, 0, 0, false, "ACGT");  // < l1, < l5
        l6.insertOccurrence(1, 1, 100, false, "ACGT");
        l6.insertOccurrence(2, 1, 100, false, "ACGT");
        l6.insertOccurrence(3, 0, 500, false, "ACGT");

        l7.insertOccurrence(0, 0, 0, false, "AGCT");  // == l1, == l5
        l7.insertOccurrence(1, 1, 100, false, "AGCT");
        l7.insertOccurrence(2, 1, 100, true, "AGCT");
        l7.insertOccurrence(3, 0, 500, false, "AGCT");

        REQUIRE(l1 == l2);
        REQUIRE(!(l1 < l2));
        REQUIRE(!(l2 < l1));
        REQUIRE(!(l1 < l3));
        REQUIRE(l3 < l1);
        REQUIRE(!(l1 == l3));
        REQUIRE(l1 < l4);
        REQUIRE(!(l4 < l1));
        REQUIRE(!(l1 == l4));
        REQUIRE(l1 == l5);
        REQUIRE(!(l5 < l1));
        REQUIRE(!(l1 < l5));
        REQUIRE(!(l1 < l6));
        REQUIRE(l6 < l1);
        REQUIRE(!(l1 == l6));
        REQUIRE(l5 == l7);
        REQUIRE(!(l7 < l5));
        REQUIRE(!(l5 < l7));
    }
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * TESTING Linkset and IdentifierMapping
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
TEST_CASE("Test Linkset on normal Testdata") {
    LoadTestdata data;
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto fastaCollection = data.fastaCollection();

    auto p = (std::thread::hardware_concurrency() > 0)
            ? std::thread::hardware_concurrency()
            : 1;
    if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }

    auto kmerMap = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort,
                                                 TwoBitKmerDataShort>>(20, genome0, genome1,
                                                                       data.idMap(),
                                                                       ULLONG_MAX, false,
                                                                       false, p);
    ExtractSeeds<TwoBitKmerDataShort,
                 TwoBitKmerDataShort>(fastaCollection,
                                      kmerMap,
                                      p, false);

    SECTION("Test with restriction max_occurrences_per_genome = 1") {
        Linkset linkset(data.idMap(), ULLONG_MAX, false, 1, 0);
        auto idMap = *(linkset.idMapping()); //constIDMapRef; // true copy so that original mapping remains unchanged

        for (auto&& elem : kmerMap->seedMap()) {
            linkset.createLinks(elem.second.at(0));
        }

        auto& lSet = linkset.linkset();

        REQUIRE(linkset.numDiscardedKmers() == 113);
        REQUIRE(linkset.numLinks() == 166);
        REQUIRE(linkset.numLinks() == lSet.size());

        // create links that are expected and check if (only) they occur
        Linkset::LinksetType expectedLinkset;
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
    auto data = LoadTestdata(LoadTestdata::geometricHashingData);
    auto genome0 = data.genome0();
    auto genome1 = data.genome1();
    auto inputFiles = data.inputFiles();
    auto fastaCollection = data.fastaCollection();

    auto p = (std::thread::hardware_concurrency() > 0)
            ? std::thread::hardware_concurrency()
            : 1;
    if (p == 1) { std::cerr << "[WARNING] -- Test only on one CPU core" << std::endl; }

    auto kmerMap = std::make_shared<SeedMapContiguous<TwoBitKmerDataShort,
                                                      TwoBitKmerDataShort>>(5, genome0, genome1,
                                                                            data.idMap(),
                                                                            ULLONG_MAX, false,
                                                                            false, p);

    ExtractSeeds<TwoBitKmerDataShort,
                 TwoBitKmerDataShort>(fastaCollection,
                                      kmerMap,
                                      p, false);

    SECTION("Test with restriction max_occurrences_per_genome = 1") {
        Linkset linkset(data.idMap(), ULLONG_MAX, false, 1, 1);
        for (auto&& elem : kmerMap->seedMap()) {
            linkset.createLinks(elem.second.at(0));
        }

        auto idMap = *(linkset.idMapping()); //constIDMapRef; // true copy so that original mapping remains unchanged
        auto& lSet = linkset.linkset();
        REQUIRE(linkset.numDiscardedKmers() == 1);
        REQUIRE(linkset.numLinks() == 15);  // as now exact positions are used in Links, TTGGC and TATGC create different Links
        REQUIRE(linkset.numLinks() == lSet.size());
        REQUIRE(idMap.refGenome() == kTestSpecies[0]);
        REQUIRE(idMap.queryGenomeID(kTestSpecies[0]) == static_cast<uint16_t>(0));
        REQUIRE(idMap.queryGenomeName(0) == kTestSpecies[0]);
        REQUIRE(idMap.querySequenceID(kTestSequences[0]) == static_cast<uint16_t>(0));
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
            auto duplicateLinkNc = std::make_shared<Link>();
            // Link 00-11 (pos. 0, 0)
            duplicateLinkNc->insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]),
                                              idMap.querySequenceID(kTestSequences[0]),
                                              0, false, "ACGT");
            duplicateLinkNc->insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]),
                                              idMap.querySequenceID(kTestSequences[1]),
                                              0, false, "ACGT");
            std::shared_ptr<Link const> duplicateLink = duplicateLinkNc;    // comparisons overloaded with <Link const>
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
        Linkset linkset(data.idMap(), ULLONG_MAX, false, ULLONG_MAX, 2);
        for (auto&& elem : kmerMap->seedMap()) {
            linkset.createLinks(elem.second.at(0));
        }

        auto idMap = *(linkset.idMapping());//constIDMapRef; // true copy so that original mapping remains unchanged
        auto& lSet = linkset.linkset();

        REQUIRE(linkset.numDiscardedKmers() == 15); // 00 - 10 from two different k-mers
        REQUIRE(linkset.numLinks() == 4);
        for (auto&& linkToCount : lSet) {
            auto count = linkToCount.second;
            REQUIRE(count == static_cast<size_t>(1));
        }

        // create links that are expected and check if (only) they occur
        auto l1 = std::make_shared<Link>();
        auto l2 = std::make_shared<Link>();
        auto l3 = std::make_shared<Link>();
        auto l4 = std::make_shared<Link>();
        // exact positions stored
        l1->insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]), idMap.querySequenceID(kTestSequences[0]), 90, false, "ACGT");   l1->insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]), idMap.querySequenceID(kTestSequences[1]), 42, false, "ACGT");
        l2->insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]), idMap.querySequenceID(kTestSequences[0]), 105, false, "ACGT");   l2->insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]), idMap.querySequenceID(kTestSequences[1]), 105, false, "ACGT");
        l3->insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]), idMap.querySequenceID(kTestSequences[0]), 90, false, "ACGT");   l3->insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]), idMap.querySequenceID(kTestSequences[1]), 105, false, "ACGT");
        l4->insertOccurrence(idMap.queryGenomeID(kTestSpecies[0]), idMap.querySequenceID(kTestSequences[0]), 105, false, "ACGT");   l4->insertOccurrence(idMap.queryGenomeID(kTestSpecies[1]), idMap.querySequenceID(kTestSequences[1]), 42, false, "ACGT");
        std::unordered_set<std::shared_ptr<Link const>,
                           LinkPtrHash, LinkPtrEqual> expectedLinkset{l1, l2, l3, l4};
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
