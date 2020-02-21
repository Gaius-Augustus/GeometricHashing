#include <algorithm>
//#include <filesystem>
//#include <fstream>
#include <iostream>

#include "catch2/catch.hpp"
#include "../FastaRepresentation.h"
#include "LoadTestdata.h"
#include "WriteTestdata.h"

TEST_CASE("FastaSequence") {
    auto seq = FastaSequence("ACGTGAGA", "sequence", "genome");
    REQUIRE(seq.genomeName() == "genome");
    REQUIRE(seq.sequenceName() == "sequence");
    REQUIRE(seq.sequence() == "ACGTGAGA");
    seq.append("CGTA");
    REQUIRE(seq.sequence() == "ACGTGAGACGTA");

    auto seq2 = FastaSequence("ACGTGAGACGTA", "sequence", "genome");
    REQUIRE(seq == seq2);

    auto seq3 = FastaSequence("ACGTGAGA", "sequence", "genome");
    REQUIRE(!(seq == seq3));
}

TEST_CASE("FastaRepresentation") {
    auto inputFileContent = ">chr1|209432133|110|+|248956422|110\n"
                            "AAGATCCTCAGACAATCCATGTGCTTCTCTTGTCCTTCATTCCACCGGAGTCTGTCTCATACCCAACCAGATTTCAGTGGAGTGAAGTTCAGGAGGCATGGAGCTGACAA\n"
                            "\n"
                            ">chr1|172138798|110|-|248956422|110\n"
                            "GGGCCTGGCTGGACAGAGTTGTCATGTGTCTGCCTGTCTACACTTGCTGTGCAGAACATCCGCTCACCTGTACAGCAGGCACAGACAGGCAGTCACATGACAACCCAGCC\n"
                            "\n";

    // create fasta file with above content
    auto testfile = WriteTestdata(inputFileContent);
    auto & tempFilename = testfile.filename();

    // create object from tmpFile
    FastaRepresentation fasta(tempFilename);
    // create object with additional artificial sequence
    FastaRepresentation fastaWithArt(tempFilename, 20);
    FastaRepresentation fastaWithArt2(tempFilename, 20); // art. sequences should differ, may fail however
    // create object with dynamic artificial sequences
    FastaRepresentation fastaWithDyn(tempFilename, 1, FastaRepresentation::dynamicallyGenerateArtificialSequences);
    FastaRepresentation fastaWithDyn2(tempFilename, 1, FastaRepresentation::dynamicallyGenerateArtificialSequences); // art. sequences should differ, may fail however

    // remove tmpFile
    testfile.deleteFile();

    REQUIRE(FastaRepresentation::genomeFromFilename("tmp_testFasta_seedFindingCpp") == "tmp_testFasta_seedFindingCpp");
    REQUIRE(fasta.genome() == tempFilename);
    REQUIRE(fastaWithArt.genome() == tempFilename);

    SECTION("Test FastaRepresentation") {
        //REQUIRE(fasta.genome() == tempFile.filename());
        REQUIRE(fasta.genome() == tempFilename);
        REQUIRE(fasta.numSequences() == 2);        
        REQUIRE(fasta.sequence("chr1|209432133|110|+|248956422|110")
                == "AAGATCCTCAGACAATCCATGTGCTTCTCTTGTCCTTCATTCCACCGGAGTCTGTCTCATACCCAACCAGATTTCAGTGGAGTGAAGTTCAGGAGGCATGGAGCTGACAA");
        REQUIRE(fasta.sequence("chr1|172138798|110|-|248956422|110")
                == "GGGCCTGGCTGGACAGAGTTGTCATGTGTCTGCCTGTCTACACTTGCTGTGCAGAACATCCGCTCACCTGTACAGCAGGCACAGACAGGCAGTCACATGACAACCCAGCC");
        REQUIRE(fasta.sequence("chr1|172138798|110|-|248956422|110")
                == fasta.fastaSequence("chr1|172138798|110|-|248956422|110").sequence());

        auto fastaHeadersObs = fasta.headers();
        auto fastaHeadersExp = std::vector<std::string>{"chr1|209432133|110|+|248956422|110",
                                                        "chr1|172138798|110|-|248956422|110"};
        for (auto&& headObs : fastaHeadersObs) {
            REQUIRE(std::find(fastaHeadersExp.begin(), fastaHeadersExp.end(), headObs) != fastaHeadersExp.end());
        }
        for (auto&& headExp : fastaHeadersExp) {
            REQUIRE(std::find(fastaHeadersObs.begin(), fastaHeadersObs.end(), headExp) != fastaHeadersObs.end());
        }

        FastaRepresentation::FastaMapType headToSeqExp;
        headToSeqExp.insert({"chr1|209432133|110|+|248956422|110", FastaSequence("AAGATCCTCAGACAATCCATGTGCTTCTCTTGTCCTTCATTCCACCGGAGTCTGTCTCATACCCAACCAGATTTCAGTGGAGTGAAGTTCAGGAGGCATGGAGCTGACAA",
                                                                                 "chr1|209432133|110|+|248956422|110", tempFilename)});
        headToSeqExp.insert({"chr1|172138798|110|-|248956422|110", FastaSequence("GGGCCTGGCTGGACAGAGTTGTCATGTGTCTGCCTGTCTACACTTGCTGTGCAGAACATCCGCTCACCTGTACAGCAGGCACAGACAGGCAGTCACATGACAACCCAGCC",
                                                                                 "chr1|172138798|110|-|248956422|110", tempFilename)});
        //headToSeqExp["chr1|172138798|110|-|248956422|110"] = "GGGCCTGGCTGGACAGAGTTGTCATGTGTCTGCCTGTCTACACTTGCTGTGCAGAACATCCGCTCACCTGTACAGCAGGCACAGACAGGCAGTCACATGACAACCCAGCC";
        REQUIRE(fasta.headerToSequence() == headToSeqExp);
   }

    SECTION("Test FastaRepresentation with artificial") {
        //REQUIRE(fastaWithArt.genome() == tempFile.filename());
        REQUIRE(fastaWithArt.genome() == tempFilename);
        REQUIRE(fastaWithArt.numSequences() == 3);
        REQUIRE(fastaWithArt.sequence("chr1|209432133|110|+|248956422|110")
                == "AAGATCCTCAGACAATCCATGTGCTTCTCTTGTCCTTCATTCCACCGGAGTCTGTCTCATACCCAACCAGATTTCAGTGGAGTGAAGTTCAGGAGGCATGGAGCTGACAA");
        REQUIRE(fastaWithArt.sequence("chr1|172138798|110|-|248956422|110")
                == "GGGCCTGGCTGGACAGAGTTGTCATGTGTCTGCCTGTCTACACTTGCTGTGCAGAACATCCGCTCACCTGTACAGCAGGCACAGACAGGCAGTCACATGACAACCCAGCC");
        REQUIRE(fastaWithArt.sequence("chr1|172138798|110|-|248956422|110")
                == fastaWithArt.fastaSequence("chr1|172138798|110|-|248956422|110").sequence());

        auto fastaWithArtHeadersObs = fastaWithArt.headers();
        auto fastaWithArtHeadersExp = std::vector<std::string>{"chr1|209432133|110|+|248956422|110",
                                                               "chr1|172138798|110|-|248956422|110",
                                                               tempFilename + "|20|artificial|0"};
        for (auto&& headObs : fastaWithArtHeadersObs) {
            REQUIRE(std::find(fastaWithArtHeadersExp.begin(), fastaWithArtHeadersExp.end(), headObs) != fastaWithArtHeadersExp.end());
        }
        for (auto&& headExp : fastaWithArtHeadersExp) {
            REQUIRE(std::find(fastaWithArtHeadersObs.begin(), fastaWithArtHeadersObs.end(), headExp) != fastaWithArtHeadersObs.end());
        }

        FastaRepresentation::FastaMapType headToSeqExp;
        headToSeqExp.insert({"chr1|209432133|110|+|248956422|110", FastaSequence("AAGATCCTCAGACAATCCATGTGCTTCTCTTGTCCTTCATTCCACCGGAGTCTGTCTCATACCCAACCAGATTTCAGTGGAGTGAAGTTCAGGAGGCATGGAGCTGACAA",
                                                                                 "chr1|209432133|110|+|248956422|110", tempFilename)});
        headToSeqExp.insert({"chr1|172138798|110|-|248956422|110", FastaSequence("GGGCCTGGCTGGACAGAGTTGTCATGTGTCTGCCTGTCTACACTTGCTGTGCAGAACATCCGCTCACCTGTACAGCAGGCACAGACAGGCAGTCACATGACAACCCAGCC",
                                                                                 "chr1|172138798|110|-|248956422|110", tempFilename)});
        headToSeqExp.insert({tempFilename + "|20|artificial|0", FastaSequence(fastaWithArt.sequence(tempFilename + "|20|artificial|0"),   // a little snakeoil here...
                                                                            tempFilename + "|20|artificial|0", tempFilename)});
        REQUIRE(fastaWithArt.headerToSequence() == headToSeqExp);

        auto artificialSeq = fastaWithArt.sequence(tempFilename + "|20|artificial|0");
        REQUIRE(artificialSeq.size() == 20);
        for (auto&& base : artificialSeq) {
            REQUIRE((base == 'A' || base == 'C' || base == 'G' || base == 'T'));
        }
        REQUIRE(artificialSeq != fastaWithArt2.sequence(tempFilename + "|20|artificial|0")); // may fail
    }

    SECTION("Test FastaRepresentation with dynamic artificial") {
        REQUIRE(fastaWithDyn.genome() == tempFilename);
        REQUIRE(fastaWithDyn.numSequences() == 4);
        REQUIRE(fastaWithDyn.sequence("chr1|209432133|110|+|248956422|110")
                == "AAGATCCTCAGACAATCCATGTGCTTCTCTTGTCCTTCATTCCACCGGAGTCTGTCTCATACCCAACCAGATTTCAGTGGAGTGAAGTTCAGGAGGCATGGAGCTGACAA");
        REQUIRE(fastaWithDyn.sequence("chr1|172138798|110|-|248956422|110")
                == "GGGCCTGGCTGGACAGAGTTGTCATGTGTCTGCCTGTCTACACTTGCTGTGCAGAACATCCGCTCACCTGTACAGCAGGCACAGACAGGCAGTCACATGACAACCCAGCC");
        REQUIRE(fastaWithDyn.sequence("chr1|172138798|110|-|248956422|110")
                == fastaWithDyn.fastaSequence("chr1|172138798|110|-|248956422|110").sequence());

        auto fastaWithDynHeadersObs = fastaWithDyn.headers();
        auto fastaWithDynHeadersExp = std::vector<std::string>{"chr1|209432133|110|+|248956422|110",
                                                               "chr1|172138798|110|-|248956422|110",
                                                               tempFilename + "|110|artificial|0",
                                                               tempFilename + "|110|artificial|1"};
        for (auto&& headObs : fastaWithDynHeadersObs) {
            REQUIRE(std::find(fastaWithDynHeadersExp.begin(), fastaWithDynHeadersExp.end(), headObs) != fastaWithDynHeadersExp.end());
        }
        for (auto&& headExp : fastaWithDynHeadersExp) {
            REQUIRE(std::find(fastaWithDynHeadersObs.begin(), fastaWithDynHeadersObs.end(), headExp) != fastaWithDynHeadersObs.end());
        }

        FastaRepresentation::FastaMapType headToSeqExp;
        headToSeqExp.insert({"chr1|209432133|110|+|248956422|110", FastaSequence("AAGATCCTCAGACAATCCATGTGCTTCTCTTGTCCTTCATTCCACCGGAGTCTGTCTCATACCCAACCAGATTTCAGTGGAGTGAAGTTCAGGAGGCATGGAGCTGACAA",
                                                                                 "chr1|209432133|110|+|248956422|110", tempFilename)});
        headToSeqExp.insert({"chr1|172138798|110|-|248956422|110", FastaSequence("GGGCCTGGCTGGACAGAGTTGTCATGTGTCTGCCTGTCTACACTTGCTGTGCAGAACATCCGCTCACCTGTACAGCAGGCACAGACAGGCAGTCACATGACAACCCAGCC",
                                                                                 "chr1|172138798|110|-|248956422|110", tempFilename)});
        headToSeqExp.insert({tempFilename + "|110|artificial|0", FastaSequence(fastaWithDyn.sequence(tempFilename + "|110|artificial|0"),   // a little snakeoil here...
                                                                               tempFilename + "|110|artificial|0", tempFilename)});
        headToSeqExp.insert({tempFilename + "|110|artificial|1", FastaSequence(fastaWithDyn.sequence(tempFilename + "|110|artificial|1"),   // a little snakeoil here...
                                                                               tempFilename + "|110|artificial|1", tempFilename)});
        REQUIRE(fastaWithDyn.headerToSequence() == headToSeqExp);

        auto artificialSeq0 = fastaWithDyn.sequence(tempFilename + "|110|artificial|0");
        auto artificialSeq1 = fastaWithDyn.sequence(tempFilename + "|110|artificial|1");
        REQUIRE(artificialSeq0.size() == 110);
        REQUIRE(artificialSeq1.size() == 110);
        for (auto&& base : artificialSeq0) {
            REQUIRE((base == 'A' || base == 'C' || base == 'G' || base == 'T'));
        }
        for (auto&& base : artificialSeq1) {
            REQUIRE((base == 'A' || base == 'C' || base == 'G' || base == 'T'));
        }
        REQUIRE(artificialSeq0 != fastaWithDyn2.sequence(tempFilename + "|110|artificial|0")); // may fail
        REQUIRE(artificialSeq1 != fastaWithDyn2.sequence(tempFilename + "|110|artificial|1")); // may fail
        REQUIRE(artificialSeq0 != artificialSeq1);   // may fail very(!) rarely
    }
}



TEST_CASE("Populate idMap") {
    auto data = LoadTestdata();
    auto idMap = data.idMap();
    auto fastaCollection = data.fastaCollection();

    REQUIRE(idMap->refGenome() == data.genome0());
    REQUIRE(idMap->queryGenomeName(0) == data.genome0());
    REQUIRE(idMap->queryGenomeName(1) == data.genome1());
    REQUIRE(idMap->queryGenomeIDConst(data.genome0()) == 0);
    REQUIRE(idMap->queryGenomeIDConst(data.genome1()) == 1);
    // not supposed to throw vvv
    auto macFasID = idMap->queryGenomeIDConst("macFas5_orthologs");
    auto hetGlaID = idMap->queryGenomeIDConst("hetGla2_orthologs");
    REQUIRE(idMap->queryGenomeName(macFasID) == "macFas5_orthologs");
    REQUIRE(idMap->queryGenomeName(hetGlaID) == "hetGla2_orthologs");
    auto seqID0 = idMap->querySequenceID("hetGla2|JH602120|2525816|67|-|10179728|67");
    auto seqID1 = idMap->querySequenceID("hetGla2|JH602084|9538844|110|+|20331017|110");
    auto seqID2 = idMap->querySequenceID("hg38|chr1|209432133|110|+|248956422|110");
    auto seqID3 = idMap->querySequenceID("hg38|chr1|172138798|110|-|248956422|110");
    auto seqID4 = idMap->querySequenceID("macFas5|chr1|30111737|110|-|227556264|110");
    auto seqID5 = idMap->querySequenceID("macFas5|chr1|67985998|109|+|227556264|109");
    auto seqID6 = idMap->querySequenceID("mm10|chr1|193507463|68|-|195471971|68");
    auto seqID7 = idMap->querySequenceID("mm10|chr1|162223368|110|+|195471971|110");
    REQUIRE(idMap->querySequenceName(seqID0) == "hetGla2|JH602120|2525816|67|-|10179728|67");
    REQUIRE(idMap->querySequenceName(seqID1) == "hetGla2|JH602084|9538844|110|+|20331017|110");
    REQUIRE(idMap->querySequenceName(seqID2) == "hg38|chr1|209432133|110|+|248956422|110");
    REQUIRE(idMap->querySequenceName(seqID3) == "hg38|chr1|172138798|110|-|248956422|110");
    REQUIRE(idMap->querySequenceName(seqID4) == "macFas5|chr1|30111737|110|-|227556264|110");
    REQUIRE(idMap->querySequenceName(seqID5) == "macFas5|chr1|67985998|109|+|227556264|109");
    REQUIRE(idMap->querySequenceName(seqID6) == "mm10|chr1|193507463|68|-|195471971|68");
    REQUIRE(idMap->querySequenceName(seqID7) == "mm10|chr1|162223368|110|+|195471971|110");
    REQUIRE(idMap->numGenomes() == 4);
    REQUIRE(idMap->numSequences() == 8);
}
