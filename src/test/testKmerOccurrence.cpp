#include <cmath>

#include "catch2/catch.hpp"
#include "../KmerOccurrence.h"

TEST_CASE("KmerOccurrence") {
    SECTION("One") {
        KmerOccurrence occ(10, 174762, 733007751850, true, "AGGTCTAGACCA");    // rc: TGGTCTAGACCT
        REQUIRE(occ.genome() == 10);
        REQUIRE(occ.sequence() == 174762);
        REQUIRE(occ.position() == 733007751850);
        REQUIRE(occ.reverse() == true);
        REQUIRE(occ.storedKmer("AGGTCTAGACCA") == "AGGTCTAGACCA");
        REQUIRE(occ.storedKmer("TGGTCTAGACCT") == "AGGTCTAGACCA");
    }
    SECTION("Two") {
        KmerOccurrence occ(5, 87381, 366503875925, false, "TGGTCTAGACCT");
        REQUIRE(occ.genome() == 5);
        REQUIRE(occ.sequence() == 87381);
        REQUIRE(occ.position() == 366503875925);
        REQUIRE(occ.reverse() == false);
        REQUIRE(occ.storedKmer("AGGTCTAGACCA") == "TGGTCTAGACCT");
        REQUIRE(occ.storedKmer("TGGTCTAGACCT") == "TGGTCTAGACCT");
    }
    SECTION("Test Genome Throw") {
        auto thrown = false;
        try {
            KmerOccurrence occ(16, 0, 0, false, "ACGT");
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);
    }
    SECTION("Test Sequence Throw") {
        auto thrown = false;
        try {
            KmerOccurrence occ(15, std::exp2(18), 0, false, "ACGT");
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);
    }
    SECTION("Test Position Throw") {
        auto thrown = false;
        try {
            KmerOccurrence occ(15, (std::exp2(18) - 1), std::exp2(40), false, "ACGT");
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);
    }
    SECTION("Test Max Position") {
        KmerOccurrence occ(0, 0, 0xffffffffff, false, "ACGT");
        REQUIRE(occ.position() == 0xffffffffff);
    }
    SECTION("Test kmer retrieval") {
        auto fastaCollection = FastaCollection();
        fastaCollection.emplace("hg38_orthologs", "./testdata/hg38_orthologs.fa");
        //auto fasta = FastaRepresentation("./testdata/hg38_orthologs.fa");
        auto & fasta = fastaCollection.fastaRepresentation("hg38_orthologs");
        auto idMap = IdentifierMapping("hg38_orthologs");
        //for (auto&& elem : fasta.headerToSequence()) {
        //    idMap.querySequenceID(elem.first, elem.second.genomeName());
        //}
        fastaCollection.populateIdentifierMappingFromFastaCollection(idMap);
        auto& seq = fasta.headerToSequence().at(idMap.querySequenceName(0)).sequence();    // first sequence
        std::unordered_map<KmerOccurrence, std::string, KmerOccurrenceHash> occToString;
        size_t k = 5;
        for (size_t i = 0; i <= (seq.size() - k); ++i) {
            auto kmer = seq.substr(i,k);
            auto kmerRC = reverseComplement(kmer);
            auto occ = KmerOccurrence(0,0,i,false,kmer);
            auto occRC = KmerOccurrence(0,0,i,true,kmerRC);
            occToString[occ] = kmer;
            occToString[occRC] = kmerRC;
        }
        for (auto&& elem : occToString) {
            auto occ = elem.first;
            auto str = elem.second;
            REQUIRE(occ.getKmerString(fastaCollection, idMap, k) == str);
        }
    }
    SECTION("Test equalSpot") {
        KmerOccurrence occ1(0,0,0,0,"TCGT");    // rc: ACGA -> c bit set
        KmerOccurrence occ2(0,0,0,0,"ACGA");    // rc: TCGT -> c bit not set
        REQUIRE(!(occ1 == occ2));
        REQUIRE(KmerOccurrence::equalSpot(occ1, occ2));
    }
    SECTION("Test centerPosition") {
        KmerOccurrence occ1(0,0,0,false,"ACGT");
        REQUIRE(KmerOccurrence::centerPosition(occ1.position(), 4) == 1);
        KmerOccurrence occ2(0,0,0,false,"ACGTA");
        REQUIRE(KmerOccurrence::centerPosition(occ2.position(), 5) == 2);
        KmerOccurrence occ3(0,0,0,false,"A");
        REQUIRE(KmerOccurrence::centerPosition(occ3.position(), 1) == 0);
    }
}

TEST_CASE("KmerOccurrencePair") {
    KmerOccurrence occ1_1(0,1,100,false,"AAAA");
    KmerOccurrence occ1_2(1,2,50,false,"AAAA");
    KmerOccurrencePair dist1(occ1_2, occ1_1);   // should store occ1 as first
    REQUIRE(dist1.first() == occ1_1);
    REQUIRE(dist1.second() == occ1_2);
    REQUIRE(dist1.distance() == -50);

    KmerOccurrence occ2_1(0,0,100,false,"AAAA");  // smaller seq
    KmerOccurrence occ2_2(1,2,50,false,"AAAA");
    KmerOccurrencePair dist2(occ2_1, occ2_2);
    REQUIRE(dist2.first() == occ2_1);
    REQUIRE(dist2.second() == occ2_2);
    REQUIRE(dist2.distance() == -50);
    REQUIRE(dist2.first() < dist1.first());
    REQUIRE(dist2.second() == dist1.second());
    REQUIRE(dist2 < dist1);
    REQUIRE(!(dist1 < dist2));
    REQUIRE(!(dist1.sameDistance(dist2)));

    KmerOccurrence occ3_1(0,1,90,false,"AAAA");  // bigger distance, smaller pos
    KmerOccurrence occ3_2(1,2,50,false,"AAAA");
    KmerOccurrencePair dist3(occ3_1, occ3_2);
    REQUIRE(dist3.distance() == -40);
    REQUIRE(dist3.first() < dist1.first());
    REQUIRE(dist3.second() == dist1.second());
    REQUIRE(dist1 < dist3); // because distance gets evaluated before positions!
    REQUIRE(!(dist1.sameDistance(dist3)));

    KmerOccurrence occ4_1(0,1,100,false,"AAAA");
    KmerOccurrence occ4_2(1,2,40,false,"AAAA");   // smaller distance, smaller pos
    KmerOccurrencePair dist4(occ4_1, occ4_2);
    REQUIRE(dist4.distance() == -60);
    REQUIRE(dist4.first() == dist1.first());
    REQUIRE(dist4.second() < dist1.second());
    REQUIRE(dist4 < dist1);
    REQUIRE(!(dist1.sameDistance(dist4)));

    KmerOccurrence occ5_1(0,1,90,false,"AAAA");  // same distance, smaller pos
    KmerOccurrence occ5_2(1,2,40,false,"AAAA");
    KmerOccurrencePair dist5(occ5_1, occ5_2);
    REQUIRE(dist5.first() < dist1.first());
    REQUIRE(dist5.second() < dist1.second());
    REQUIRE(dist1.sameDistance(dist5));
    REQUIRE(dist5.sameDistance(dist1));

    SECTION("Test Diagonal Sorting") {
        KmerOccurrencePair diag1_1(KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{1,1,10,false,"AAAA"});
        KmerOccurrencePair diag1_2(KmerOccurrence{0,0,1,false,"AAAA"}, KmerOccurrence{1,1,11,false,"AAAA"});
        KmerOccurrencePair diag1_3(KmerOccurrence{0,0,2,false,"AAAA"}, KmerOccurrence{1,1,12,false,"AAAA"});
        KmerOccurrencePair diag2_1(KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{1,1,0,false,"AAAA"});
        KmerOccurrencePair diag2_2(KmerOccurrence{0,0,1,false,"AAAA"}, KmerOccurrence{1,1,1,false,"AAAA"});
        KmerOccurrencePair diag2_3(KmerOccurrence{0,0,2,false,"AAAA"}, KmerOccurrence{1,1,2,false,"AAAA"});
        KmerOccurrencePair diag3_1(KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{2,2,10,false,"AAAA"});
        KmerOccurrencePair diag3_2(KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{2,2,11,false,"AAAA"});
        KmerOccurrencePair diag3_3(KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{2,2,12,false,"AAAA"});
        std::vector<KmerOccurrencePair> unsorted{diag3_3, diag2_1, diag1_2,
                                                     diag1_1, diag2_3, diag1_3,
                                                     diag2_2, diag3_2, diag3_1};
        std::vector<KmerOccurrencePair> sorted{diag2_1, diag2_2, diag2_3,
                                                   diag1_1, diag1_2, diag1_3,
                                                   diag3_1, diag3_2, diag3_3};
        std::sort(unsorted.begin(), unsorted.end());
        REQUIRE(unsorted == sorted);
    }
}
