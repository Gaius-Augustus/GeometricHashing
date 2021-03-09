#include <cmath>

#include "catch2/catch.hpp"
#include "../KmerOccurrence.h"
#include "LoadTestdata.h"

TEST_CASE("KmerOccurrence") {
    std::cout << "[INFO] -- [TEST CASE] -- KmerOccurrence" << std::endl;
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
        auto data = LoadTestdata();
        fastaCollection.emplace("hg38_orthologs", data.inputFiles().at(1));
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
