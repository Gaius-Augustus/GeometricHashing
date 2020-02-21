#include "catch2/catch.hpp"

#include "../IdentifierMapping.h"

TEST_CASE("Identifier Mapping") {
    IdentifierMapping idMap("referenceGenome");
    // test reference creation
    REQUIRE(idMap.queryGenomeIDConst("referenceGenome") == 0);
    REQUIRE(idMap.queryGenomeID("referenceGenome") == 0);
    REQUIRE(idMap.queryGenomeName(0) == "referenceGenome");
    REQUIRE(idMap.refGenome() == "referenceGenome");
    REQUIRE(idMap.numGenomes() == 1);
    REQUIRE(idMap.numSequences() == 0);
    //REQUIRE(idMap.referenceGenomeIDQueried() == 2);
    // test genome insertion
    idMap.queryGenomeID("secondGenome");
    REQUIRE(idMap.queryGenomeIDConst("secondGenome") == 1);
    REQUIRE(idMap.queryGenomeID("secondGenome") == 1);
    REQUIRE(idMap.queryGenomeName(1) == "secondGenome");
    REQUIRE(idMap.numGenomes() == 2);
    REQUIRE(idMap.numSequences() == 0);
    //REQUIRE(idMap.referenceGenomeIDQueried() == 2);
    // test sequence insertion
    idMap.querySequenceID("sequence0", "referenceGenome");  // referenceGenomeIDQueried += 1
    idMap.querySequenceID("sequence1", "secondGenome");
    idMap.querySequenceID("sequence2", "referenceGenome");  // referenceGenomeIDQueried += 1
    REQUIRE(idMap.querySequenceID("sequence0") == 0);
    REQUIRE(idMap.querySequenceID("sequence0", "referenceGenome") == 0);
    REQUIRE(idMap.querySequenceName(0) == "sequence0");
    REQUIRE(idMap.sequenceIDToGenomeID().at(0) == 0);
    REQUIRE(idMap.querySequenceID("sequence1") == 1);
    REQUIRE(idMap.querySequenceName(1) == "sequence1");
    REQUIRE(idMap.sequenceIDToGenomeID().at(1) == 1);
    REQUIRE(idMap.querySequenceID("sequence2") == 2);
    REQUIRE(idMap.querySequenceName(2) == "sequence2");
    REQUIRE(idMap.sequenceIDToGenomeID().at(2) == 0);
    REQUIRE(idMap.numGenomes() == 2);
    REQUIRE(idMap.numSequences() == 3);
    //REQUIRE(idMap.referenceGenomeIDQueried() == 4);
    // test throwing
    auto thrown = false;
    try {
        idMap.queryGenomeIDConst("thirdGenome");
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        idMap.queryGenomeName(2);
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        idMap.querySequenceID("sequence3");
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        idMap.querySequenceName(3);
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
    REQUIRE(idMap2.refGenome() == "referenceGenome");
    REQUIRE(idMap2.numGenomes() == 2);
    REQUIRE(idMap2.numSequences() == 3);
    //REQUIRE(idMap2.referenceGenomeIDQueried() == 4);
    REQUIRE(idMap2.queryGenomeIDConst("referenceGenome") == 0); // referenceGenomeIDQueried += 1
    REQUIRE(idMap2.queryGenomeIDConst("secondGenome") == 1);
    REQUIRE(idMap2.querySequenceID("sequence0") == 0);
    REQUIRE(idMap2.querySequenceID("sequence1") == 1);
    REQUIRE(idMap2.querySequenceID("sequence2") == 2);
    //REQUIRE(idMap2.referenceGenomeIDQueried() == 5);
    // check if everything remained
    //REQUIRE(idMap.referenceGenomeIDQueried() == 4);
    REQUIRE(idMap.queryGenomeIDConst("referenceGenome") == 0);
    REQUIRE(idMap.queryGenomeName(0) == "referenceGenome");
    REQUIRE(idMap.refGenome() == "referenceGenome");
    REQUIRE(idMap.numGenomes() == 2);
    REQUIRE(idMap.numSequences() == 3);
    //REQUIRE(idMap.referenceGenomeIDQueried() == 5);
}
