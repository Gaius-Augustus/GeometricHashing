#include "catch2/catch.hpp"

#include "../IdentifierMapping.h"

TEST_CASE("Identifier Mapping") {
    std::cout << "[INFO] -- [TEST CASE] -- Identifier Mapping" << std::endl;
    IdentifierMapping idMap("referenceGenome");
    // test reference creation
    REQUIRE(idMap.queryGenomeIDConst("referenceGenome") == 0);
    REQUIRE(idMap.queryGenomeID("referenceGenome") == 0);
    REQUIRE(idMap.queryGenomeName(0) == "referenceGenome");
    //REQUIRE(idMap.refGenome() == "referenceGenome");
    REQUIRE(idMap.numGenomes() == 1);
    REQUIRE(idMap.numSequences() == 0);
    // test if counters are correct
    REQUIRE(idMap.genomeIDToName().size() == idMap.genomeNameToID().size());
    REQUIRE(idMap.genomeIDToName().size() == idMap.numGenomes());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.sequenceTupleToID().size());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.numSequences());
    //REQUIRE(idMap.referenceGenomeIDQueried() == 2);
    // test genome insertion
    idMap.queryGenomeID("secondGenome");
    REQUIRE(idMap.queryGenomeIDConst("secondGenome") == 1);
    REQUIRE(idMap.queryGenomeID("secondGenome") == 1);
    REQUIRE(idMap.queryGenomeName(1) == "secondGenome");
    REQUIRE(idMap.numGenomes() == 2);
    REQUIRE(idMap.numSequences() == 0);
    // test if counters are correct
    REQUIRE(idMap.genomeIDToName().size() == idMap.genomeNameToID().size());
    REQUIRE(idMap.genomeIDToName().size() == idMap.numGenomes());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.sequenceTupleToID().size());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.numSequences());
    //REQUIRE(idMap.referenceGenomeIDQueried() == 2);
    // test sequence insertion
    idMap.querySequenceID("sequence0", "referenceGenome");  // referenceGenomeIDQueried += 1
    // test if counters are correct
    REQUIRE(idMap.genomeIDToName().size() == idMap.genomeNameToID().size());
    REQUIRE(idMap.genomeIDToName().size() == idMap.numGenomes());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.sequenceTupleToID().size());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.numSequences());
    idMap.querySequenceID("sequence1", "secondGenome");
    // test if counters are correct
    REQUIRE(idMap.genomeIDToName().size() == idMap.genomeNameToID().size());
    REQUIRE(idMap.genomeIDToName().size() == idMap.numGenomes());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.sequenceTupleToID().size());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.numSequences());
    idMap.querySequenceID("sequence2", "referenceGenome");  // referenceGenomeIDQueried += 1
    // test if counters are correct
    REQUIRE(idMap.genomeIDToName().size() == idMap.genomeNameToID().size());
    REQUIRE(idMap.genomeIDToName().size() == idMap.numGenomes());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.sequenceTupleToID().size());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.numSequences());
    idMap.querySequenceID("sequence1", "thirdGenome");
    // test if counters are correct
    REQUIRE(idMap.genomeIDToName().size() == idMap.genomeNameToID().size());
    REQUIRE(idMap.genomeIDToName().size() == idMap.numGenomes());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.sequenceTupleToID().size());
    REQUIRE(idMap.sequenceIDToTuple().size() == idMap.numSequences());

    REQUIRE(idMap.querySequenceIDConst("sequence0", "referenceGenome") == 0);
    REQUIRE(idMap.querySequenceID("sequence0", "referenceGenome") == 0);
    REQUIRE(idMap.querySequenceName(0) == "sequence0");
    REQUIRE(idMap.sequenceIDToTuple().at(0) == IdentifierMapping::SequenceGenomeTuple{0, "sequence0"});
    REQUIRE(idMap.genomeIDFromSequenceID(0) == 0);
    REQUIRE(idMap.querySequenceIDConst("sequence1", "secondGenome") == 1);
    REQUIRE(idMap.querySequenceID("sequence1", "secondGenome") == 1);
    REQUIRE(idMap.querySequenceName(1) == "sequence1");
    REQUIRE(idMap.sequenceIDToTuple().at(1) == IdentifierMapping::SequenceGenomeTuple{1, "sequence1"});
    REQUIRE(idMap.genomeIDFromSequenceID(1) == 1);
    REQUIRE(idMap.querySequenceIDConst("sequence2", "referenceGenome") == 2);
    REQUIRE(idMap.querySequenceID("sequence2", "referenceGenome") == 2);
    REQUIRE(idMap.querySequenceName(2) == "sequence2");
    REQUIRE(idMap.sequenceIDToTuple().at(2) == IdentifierMapping::SequenceGenomeTuple{0, "sequence2"});
    REQUIRE(idMap.genomeIDFromSequenceID(2) == 0);
    REQUIRE(idMap.querySequenceIDConst("sequence1", "thirdGenome") == 3);
    REQUIRE(idMap.querySequenceID("sequence1", "thirdGenome") == 3);
    REQUIRE(idMap.querySequenceName(3) == "sequence1");
    REQUIRE(idMap.sequenceIDToTuple().at(3) == IdentifierMapping::SequenceGenomeTuple{2, "sequence1"});
    REQUIRE(idMap.genomeIDFromSequenceID(3) == 2);

    REQUIRE(idMap.numGenomes() == 3);
    REQUIRE(idMap.numSequences() == 4);
    //REQUIRE(idMap.referenceGenomeIDQueried() == 4);
    // test throwing
    auto thrown = false;
    try {
        idMap.queryGenomeIDConst("fourthGenome");
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        idMap.queryGenomeName(3);
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        idMap.querySequenceIDConst("sequence3", "fourthGenome");
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        idMap.querySequenceIDConst("sequence3", "thirdGenome");
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        idMap.querySequenceIDConst("sequence1", "fourthGenome");
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        idMap.querySequenceName(4);
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    // check copying
    IdentifierMapping idMap2("something");
    idMap2.completeUpdate(idMap);
    REQUIRE(idMap == idMap2);
    // check if idMap2 is ok
    //REQUIRE(idMap2.refGenome() == "referenceGenome");
    REQUIRE(idMap2.numGenomes() == 3);
    REQUIRE(idMap2.numSequences() == 4);
    //REQUIRE(idMap2.referenceGenomeIDQueried() == 4);
    REQUIRE(idMap2.queryGenomeIDConst("referenceGenome") == 0); // referenceGenomeIDQueried += 1
    REQUIRE(idMap2.queryGenomeIDConst("secondGenome") == 1);
    REQUIRE(idMap2.queryGenomeIDConst("thirdGenome") == 2);
    REQUIRE(idMap2.querySequenceIDConst("sequence0", "referenceGenome") == 0);
    REQUIRE(idMap2.querySequenceIDConst("sequence1", "secondGenome") == 1);
    REQUIRE(idMap2.querySequenceIDConst("sequence2", "referenceGenome") == 2);
    REQUIRE(idMap2.querySequenceIDConst("sequence1", "thirdGenome") == 3);
    //REQUIRE(idMap2.referenceGenomeIDQueried() == 5);
    // check if everything remained
    //REQUIRE(idMap.referenceGenomeIDQueried() == 4);
    REQUIRE(idMap.queryGenomeIDConst("referenceGenome") == 0);
    REQUIRE(idMap.queryGenomeName(0) == "referenceGenome");
    //REQUIRE(idMap.refGenome() == "referenceGenome");
    REQUIRE(idMap.numGenomes() == 3);
    REQUIRE(idMap.numSequences() == 4);
    //REQUIRE(idMap.referenceGenomeIDQueried() == 5);
}
