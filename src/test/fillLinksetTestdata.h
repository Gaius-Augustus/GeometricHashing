#ifndef FILLLINKSETTESTDATA_H
#define FILLLINKSETTESTDATA_H

#include "../IdentifierMapping.h"
#include "../Link.h"
#include "../Linkset.h"

inline void fillLinksetTestdata_part0(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link0 = Link();
    link0.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 110, false, "TTCAGGAGGCATGGAGCTGA");
    link0.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 109, false, "TTCAGGAGGCATGGAGCTGA");
    if (linkset.find(link0) != linkset.end()) { ++linkset.at(link0); } else { linkset[link0] = 1; }
    auto link1 = Link();
    link1.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 38, false, "CACTTGCTGTGCAGAACATC");
    link1.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 40, false, "CACTTGCTGTGCAGAACATC");
    link1.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 38, false, "CACTTGCTGTGCAGAACATC");
    link1.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 40, false, "CACTTGCTGTGCAGAACATC");
    if (linkset.find(link1) != linkset.end()) { ++linkset.at(link1); } else { linkset[link1] = 1; }
    auto link2 = Link();
    link2.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 109, false, "GTTCAGGAGGCATGGAGCTG");
    link2.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 108, false, "GTTCAGGAGGCATGGAGCTG");
    if (linkset.find(link2) != linkset.end()) { ++linkset.at(link2); } else { linkset[link2] = 1; }
    auto link3 = Link();
    link3.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 20, false, "GTGCTTCTCTTGTCCTTCAT");
    link3.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 20, false, "GTGCTTCTCTTGTCCTTCAT");
    if (linkset.find(link3) != linkset.end()) { ++linkset.at(link3); } else { linkset[link3] = 1; }
    auto link4 = Link();
    link4.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 64, false, "AACCAGATTTCAGTGGAGTG");
    link4.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 38, false, "AACCAGATTTCAGTGGAGTG");
    link4.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 64, false, "AACCAGATTTCAGTGGAGTG");
    if (linkset.find(link4) != linkset.end()) { ++linkset.at(link4); } else { linkset[link4] = 1; }
    auto link5 = Link();
    link5.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 32, false, "TGTCTACACTTGCTGTGCAG");
    link5.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 34, false, "TGTCTACACTTGCTGTGCAG");
    link5.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 32, false, "TGTCTACACTTGCTGTGCAG");
    link5.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 34, false, "TGTCTACACTTGCTGTGCAG");
    if (linkset.find(link5) != linkset.end()) { ++linkset.at(link5); } else { linkset[link5] = 1; }
    auto link6 = Link();
    link6.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 71, false, "AGCAGGCACAGACAGGCAGT");
    link6.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 73, false, "AGCAGGCACAGACAGGCAGT");
    link6.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 71, false, "AGCAGGCACAGACAGGCAGT");
    link6.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 73, false, "AGCAGGCACAGACAGGCAGT");
    if (linkset.find(link6) != linkset.end()) { ++linkset.at(link6); } else { linkset[link6] = 1; }
    auto link7 = Link();
    link7.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 29, false, "TTGTCCTTCATTCCACCGGA");
    link7.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 4, false, "TTGTCCTTCATTCCACCGGA");
    link7.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 4, false, "TTGTCCTTCATTCCACCGGA");
    link7.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 29, false, "TTGTCCTTCATTCCACCGGA");
    if (linkset.find(link7) != linkset.end()) { ++linkset.at(link7); } else { linkset[link7] = 1; }
    auto link8 = Link();
    link8.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 87, false, "CAGTCACATGACAACCCAGC");
    link8.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 89, false, "CAGTCACATGACAACCCAGC");
    link8.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 87, false, "CAGTCACATGACAACCCAGC");
    link8.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 89, false, "CAGTCACATGACAACCCAGC");
    if (linkset.find(link8) != linkset.end()) { ++linkset.at(link8); } else { linkset[link8] = 1; }
    auto link9 = Link();
    link9.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 22, false, "GCTTCTCTTGTCCTTCATTC");
    link9.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 22, false, "GCTTCTCTTGTCCTTCATTC");
    if (linkset.find(link9) != linkset.end()) { ++linkset.at(link9); } else { linkset[link9] = 1; }
    auto link10 = Link();
    link10.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 92, false, "GATTTCAGTGGAGTGAAGTT");
    link10.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 69, false, "GATTTCAGTGGAGTGAAGTT");
    if (linkset.find(link10) != linkset.end()) { ++linkset.at(link10); } else { linkset[link10] = 1; }
    auto link11 = Link();
    link11.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 25, false, "GTCTGCCTGTCTACACTTGC");
    link11.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 27, false, "GTCTGCCTGTCTACACTTGC");
    link11.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 25, false, "GTCTGCCTGTCTACACTTGC");
    link11.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 27, false, "GTCTGCCTGTCTACACTTGC");
    if (linkset.find(link11) != linkset.end()) { ++linkset.at(link11); } else { linkset[link11] = 1; }
    auto link12 = Link();
    link12.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 72, false, "GCAGGCACAGACAGGCAGTC");
    link12.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 74, false, "GCAGGCACAGACAGGCAGTC");
    link12.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 72, false, "GCAGGCACAGACAGGCAGTC");
    link12.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 74, false, "GCAGGCACAGACAGGCAGTC");
    if (linkset.find(link12) != linkset.end()) { ++linkset.at(link12); } else { linkset[link12] = 1; }
    auto link13 = Link();
    link13.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 51, false, "CTGTCTCATACCCAACCAGA");
    link13.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 51, false, "CTGTCTCATACCCAACCAGA");
    if (linkset.find(link13) != linkset.end()) { ++linkset.at(link13); } else { linkset[link13] = 1; }
    auto link14 = Link();
    link14.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 50, false, "AGAACATCCGCTCACCTGTA");
    link14.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 52, false, "AGAACATCCGCTCACCTGTA");
    link14.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 50, false, "AGAACATCCGCTCACCTGTA");
    link14.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 52, false, "AGAACATCCGCTCACCTGTA");
    if (linkset.find(link14) != linkset.end()) { ++linkset.at(link14); } else { linkset[link14] = 1; }
    auto link15 = Link();
    link15.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 68, false, "AGATTTCAGTGGAGTGAAGT");
    link15.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 68, false, "AGATTTCAGTGGAGTGAAGT");
    if (linkset.find(link15) != linkset.end()) { ++linkset.at(link15); } else { linkset[link15] = 1; }
    auto link16 = Link();
    link16.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 0, false, "GCCTGGCTGGACAGAGTTGT");
    link16.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 2, false, "GCCTGGCTGGACAGAGTTGT");
    link16.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 0, false, "GCCTGGCTGGACAGAGTTGT");
    link16.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 2, false, "GCCTGGCTGGACAGAGTTGT");
    if (linkset.find(link16) != linkset.end()) { ++linkset.at(link16); } else { linkset[link16] = 1; }
    auto link17 = Link();
    link17.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 13, false, "GAGTTGTCATGTGTCTGCCT");
    link17.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 15, false, "GAGTTGTCATGTGTCTGCCT");
    link17.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 13, false, "GAGTTGTCATGTGTCTGCCT");
    link17.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 15, false, "GAGTTGTCATGTGTCTGCCT");
    if (linkset.find(link17) != linkset.end()) { ++linkset.at(link17); } else { linkset[link17] = 1; }
    auto link18 = Link();
    link18.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 26, false, "CTCTTGTCCTTCATTCCACC");
    link18.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 1, false, "CTCTTGTCCTTCATTCCACC");
    link18.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 1, false, "CTCTTGTCCTTCATTCCACC");
    link18.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 26, false, "CTCTTGTCCTTCATTCCACC");
    if (linkset.find(link18) != linkset.end()) { ++linkset.at(link18); } else { linkset[link18] = 1; }
    auto link19 = Link();
    link19.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 45, false, "CGGAGTCTGTCTCATACCCA");
    link19.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 45, false, "CGGAGTCTGTCTCATACCCA");
    if (linkset.find(link19) != linkset.end()) { ++linkset.at(link19); } else { linkset[link19] = 1; }
}

inline void fillLinksetTestdata_part1(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link20 = Link();
    link20.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 27, false, "CTGCCTGTCTACACTTGCTG");
    link20.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 29, false, "CTGCCTGTCTACACTTGCTG");
    link20.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 27, false, "CTGCCTGTCTACACTTGCTG");
    link20.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 29, false, "CTGCCTGTCTACACTTGCTG");
    if (linkset.find(link20) != linkset.end()) { ++linkset.at(link20); } else { linkset[link20] = 1; }
    auto link21 = Link();
    link21.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 18, false, "ATGTGCTTCTCTTGTCCTTC");
    link21.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 18, false, "ATGTGCTTCTCTTGTCCTTC");
    if (linkset.find(link21) != linkset.end()) { ++linkset.at(link21); } else { linkset[link21] = 1; }
    auto link22 = Link();
    link22.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 12, false, "AGAGTTGTCATGTGTCTGCC");
    link22.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 14, false, "AGAGTTGTCATGTGTCTGCC");
    link22.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 12, false, "AGAGTTGTCATGTGTCTGCC");
    link22.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 14, false, "AGAGTTGTCATGTGTCTGCC");
    if (linkset.find(link22) != linkset.end()) { ++linkset.at(link22); } else { linkset[link22] = 1; }
    auto link23 = Link();
    link23.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 28, false, "CTTGTCCTTCATTCCACCGG");
    link23.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 3, false, "CTTGTCCTTCATTCCACCGG");
    link23.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 3, false, "CTTGTCCTTCATTCCACCGG");
    link23.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 28, false, "CTTGTCCTTCATTCCACCGG");
    if (linkset.find(link23) != linkset.end()) { ++linkset.at(link23); } else { linkset[link23] = 1; }
    auto link24 = Link();
    link24.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 61, false, "CCCAACCAGATTTCAGTGGA");
    link24.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 61, false, "CCCAACCAGATTTCAGTGGA");
    if (linkset.find(link24) != linkset.end()) { ++linkset.at(link24); } else { linkset[link24] = 1; }
    auto link25 = Link();
    link25.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 1, false, "CCTGGCTGGACAGAGTTGTC");
    link25.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 3, false, "CCTGGCTGGACAGAGTTGTC");
    link25.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 1, false, "CCTGGCTGGACAGAGTTGTC");
    link25.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 3, false, "CCTGGCTGGACAGAGTTGTC");
    if (linkset.find(link25) != linkset.end()) { ++linkset.at(link25); } else { linkset[link25] = 1; }
    auto link26 = Link();
    link26.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 19, false, "TGTGCTTCTCTTGTCCTTCA");
    link26.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 19, false, "TGTGCTTCTCTTGTCCTTCA");
    if (linkset.find(link26) != linkset.end()) { ++linkset.at(link26); } else { linkset[link26] = 1; }
    auto link27 = Link();
    link27.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 10, false, "ACAGAGTTGTCATGTGTCTG");
    link27.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 12, false, "ACAGAGTTGTCATGTGTCTG");
    link27.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 10, false, "ACAGAGTTGTCATGTGTCTG");
    link27.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 12, false, "ACAGAGTTGTCATGTGTCTG");
    if (linkset.find(link27) != linkset.end()) { ++linkset.at(link27); } else { linkset[link27] = 1; }
    auto link28 = Link();
    link28.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 54, false, "CATCCGCTCACCTGTACAGC");
    link28.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 56, false, "CATCCGCTCACCTGTACAGC");
    link28.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 54, false, "CATCCGCTCACCTGTACAGC");
    link28.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 56, false, "CATCCGCTCACCTGTACAGC");
    if (linkset.find(link28) != linkset.end()) { ++linkset.at(link28); } else { linkset[link28] = 1; }
    auto link29 = Link();
    link29.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 83, false, "CAGGCAGTCACATGACAACC");
    link29.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 85, false, "CAGGCAGTCACATGACAACC");
    link29.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 83, false, "CAGGCAGTCACATGACAACC");
    link29.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 85, false, "CAGGCAGTCACATGACAACC");
    if (linkset.find(link29) != linkset.end()) { ++linkset.at(link29); } else { linkset[link29] = 1; }
    auto link30 = Link();
    link30.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 20, false, "CATGTGTCTGCCTGTCTACA");
    link30.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 22, false, "CATGTGTCTGCCTGTCTACA");
    link30.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 20, false, "CATGTGTCTGCCTGTCTACA");
    link30.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 22, false, "CATGTGTCTGCCTGTCTACA");
    if (linkset.find(link30) != linkset.end()) { ++linkset.at(link30); } else { linkset[link30] = 1; }
    auto link31 = Link();
    link31.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 27, false, "TCTTGTCCTTCATTCCACCG");
    link31.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 2, false, "TCTTGTCCTTCATTCCACCG");
    link31.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 2, false, "TCTTGTCCTTCATTCCACCG");
    link31.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 27, false, "TCTTGTCCTTCATTCCACCG");
    if (linkset.find(link31) != linkset.end()) { ++linkset.at(link31); } else { linkset[link31] = 1; }
    auto link32 = Link();
    link32.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 30, false, "CCTGTCTACACTTGCTGTGC");
    link32.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 32, false, "CCTGTCTACACTTGCTGTGC");
    link32.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 30, false, "CCTGTCTACACTTGCTGTGC");
    link32.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 32, false, "CCTGTCTACACTTGCTGTGC");
    if (linkset.find(link32) != linkset.end()) { ++linkset.at(link32); } else { linkset[link32] = 1; }
    auto link33 = Link();
    link33.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 15, false, "GTTGTCATGTGTCTGCCTGT");
    link33.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 17, false, "GTTGTCATGTGTCTGCCTGT");
    link33.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 15, false, "GTTGTCATGTGTCTGCCTGT");
    link33.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 17, false, "GTTGTCATGTGTCTGCCTGT");
    if (linkset.find(link33) != linkset.end()) { ++linkset.at(link33); } else { linkset[link33] = 1; }
    auto link34 = Link();
    link34.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 53, false, "ACATCCGCTCACCTGTACAG");
    link34.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 55, false, "ACATCCGCTCACCTGTACAG");
    link34.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 53, false, "ACATCCGCTCACCTGTACAG");
    link34.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 55, false, "ACATCCGCTCACCTGTACAG");
    if (linkset.find(link34) != linkset.end()) { ++linkset.at(link34); } else { linkset[link34] = 1; }
    auto link35 = Link();
    link35.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 41, false, "CCACCGGAGTCTGTCTCATA");
    link35.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 41, false, "CCACCGGAGTCTGTCTCATA");
    if (linkset.find(link35) != linkset.end()) { ++linkset.at(link35); } else { linkset[link35] = 1; }
    auto link36 = Link();
    link36.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 18, false, "GTCATGTGTCTGCCTGTCTA");
    link36.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 20, false, "GTCATGTGTCTGCCTGTCTA");
    link36.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 18, false, "GTCATGTGTCTGCCTGTCTA");
    link36.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 20, false, "GTCATGTGTCTGCCTGTCTA");
    if (linkset.find(link36) != linkset.end()) { ++linkset.at(link36); } else { linkset[link36] = 1; }
    auto link37 = Link();
    link37.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 59, false, "TACCCAACCAGATTTCAGTG");
    link37.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 59, false, "TACCCAACCAGATTTCAGTG");
    if (linkset.find(link37) != linkset.end()) { ++linkset.at(link37); } else { linkset[link37] = 1; }
    auto link38 = Link();
    link38.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 5, false, "GCTGGACAGAGTTGTCATGT");
    link38.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 7, false, "GCTGGACAGAGTTGTCATGT");
    link38.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 5, false, "GCTGGACAGAGTTGTCATGT");
    link38.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 7, false, "GCTGGACAGAGTTGTCATGT");
    if (linkset.find(link38) != linkset.end()) { ++linkset.at(link38); } else { linkset[link38] = 1; }
}

inline void fillLinksetTestdata_part2(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link39 = Link();
    link39.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 49, false, "GTCTGTCTCATACCCAACCA");
    link39.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 49, false, "GTCTGTCTCATACCCAACCA");
    if (linkset.find(link39) != linkset.end()) { ++linkset.at(link39); } else { linkset[link39] = 1; }
    auto link40 = Link();
    link40.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 31, false, "CTGTCTACACTTGCTGTGCA");
    link40.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 33, false, "CTGTCTACACTTGCTGTGCA");
    link40.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 31, false, "CTGTCTACACTTGCTGTGCA");
    link40.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 33, false, "CTGTCTACACTTGCTGTGCA");
    if (linkset.find(link40) != linkset.end()) { ++linkset.at(link40); } else { linkset[link40] = 1; }
    auto link41 = Link();
    link41.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 58, false, "ATACCCAACCAGATTTCAGT");
    link41.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 58, false, "ATACCCAACCAGATTTCAGT");
    if (linkset.find(link41) != linkset.end()) { ++linkset.at(link41); } else { linkset[link41] = 1; }
    auto link42 = Link();
    link42.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 47, false, "TGCAGAACATCCGCTCACCT");
    link42.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 49, false, "TGCAGAACATCCGCTCACCT");
    link42.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 47, false, "TGCAGAACATCCGCTCACCT");
    link42.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 49, false, "TGCAGAACATCCGCTCACCT");
    if (linkset.find(link42) != linkset.end()) { ++linkset.at(link42); } else { linkset[link42] = 1; }
    auto link43 = Link();
    link43.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 77, false, "CACAGACAGGCAGTCACATG");
    link43.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 79, false, "CACAGACAGGCAGTCACATG");
    link43.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 77, false, "CACAGACAGGCAGTCACATG");
    link43.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 79, false, "CACAGACAGGCAGTCACATG");
    if (linkset.find(link43) != linkset.end()) { ++linkset.at(link43); } else { linkset[link43] = 1; }
    auto link44 = Link();
    link44.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 23, false, "CTTCTCTTGTCCTTCATTCC");
    link44.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 23, false, "CTTCTCTTGTCCTTCATTCC");
    if (linkset.find(link44) != linkset.end()) { ++linkset.at(link44); } else { linkset[link44] = 1; }
    auto link45 = Link();
    link45.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 60, false, "CTCACCTGTACAGCAGGCAC");
    link45.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 62, false, "CTCACCTGTACAGCAGGCAC");
    link45.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 60, false, "CTCACCTGTACAGCAGGCAC");
    link45.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 62, false, "CTCACCTGTACAGCAGGCAC");
    if (linkset.find(link45) != linkset.end()) { ++linkset.at(link45); } else { linkset[link45] = 1; }
    auto link46 = Link();
    link46.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 39, false, "TTCCACCGGAGTCTGTCTCA");
    link46.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 39, false, "TTCCACCGGAGTCTGTCTCA");
    if (linkset.find(link46) != linkset.end()) { ++linkset.at(link46); } else { linkset[link46] = 1; }
    auto link47 = Link();
    link47.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 103, false, "AGTGAAGTTCAGGAGGCATG");
    link47.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 102, false, "AGTGAAGTTCAGGAGGCATG");
    if (linkset.find(link47) != linkset.end()) { ++linkset.at(link47); } else { linkset[link47] = 1; }
    auto link48 = Link();
    link48.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 95, false, "TTCAGTGGAGTGAAGTTCAG");
    link48.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 94, false, "TTCAGTGGAGTGAAGTTCAG");
    if (linkset.find(link48) != linkset.end()) { ++linkset.at(link48); } else { linkset[link48] = 1; }
    auto link49 = Link();
    link49.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 100, false, "TGGAGTGAAGTTCAGGAGGC");
    link49.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 99, false, "TGGAGTGAAGTTCAGGAGGC");
    if (linkset.find(link49) != linkset.end()) { ++linkset.at(link49); } else { linkset[link49] = 1; }
    auto link50 = Link();
    link50.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 40, false, "CTTGCTGTGCAGAACATCCG");
    link50.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 42, false, "CTTGCTGTGCAGAACATCCG");
    link50.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 40, false, "CTTGCTGTGCAGAACATCCG");
    link50.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 42, false, "CTTGCTGTGCAGAACATCCG");
    if (linkset.find(link50) != linkset.end()) { ++linkset.at(link50); } else { linkset[link50] = 1; }
    auto link51 = Link();
    link51.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 21, false, "TGCTTCTCTTGTCCTTCATT");
    link51.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 21, false, "TGCTTCTCTTGTCCTTCATT");
    if (linkset.find(link51) != linkset.end()) { ++linkset.at(link51); } else { linkset[link51] = 1; }
    auto link52 = Link();
    link52.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 105, false, "TGAAGTTCAGGAGGCATGGA");
    link52.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 104, false, "TGAAGTTCAGGAGGCATGGA");
    if (linkset.find(link52) != linkset.end()) { ++linkset.at(link52); } else { linkset[link52] = 1; }
    auto link53 = Link();
    link53.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 16, false, "CCATGTGCTTCTCTTGTCCT");
    link53.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 16, false, "CCATGTGCTTCTCTTGTCCT");
    if (linkset.find(link53) != linkset.end()) { ++linkset.at(link53); } else { linkset[link53] = 1; }
    auto link54 = Link();
    link54.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 82, false, "ACAGGCAGTCACATGACAAC");
    link54.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 84, false, "ACAGGCAGTCACATGACAAC");
    link54.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 82, false, "ACAGGCAGTCACATGACAAC");
    link54.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 84, false, "ACAGGCAGTCACATGACAAC");
    if (linkset.find(link54) != linkset.end()) { ++linkset.at(link54); } else { linkset[link54] = 1; }
    auto link55 = Link();
    link55.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 8, false, "GGACAGAGTTGTCATGTGTC");
    link55.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 10, false, "GGACAGAGTTGTCATGTGTC");
    link55.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 8, false, "GGACAGAGTTGTCATGTGTC");
    link55.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 10, false, "GGACAGAGTTGTCATGTGTC");
    if (linkset.find(link55) != linkset.end()) { ++linkset.at(link55); } else { linkset[link55] = 1; }
    auto link56 = Link();
    link56.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 111, false, "TCAGGAGGCATGGAGCTGAC");
    link56.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 110, false, "TCAGGAGGCATGGAGCTGAC");
    if (linkset.find(link56) != linkset.end()) { ++linkset.at(link56); } else { linkset[link56] = 1; }
    auto link57 = Link();
    link57.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 36, false, "TCATTCCACCGGAGTCTGTC");
    link57.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 11, false, "TCATTCCACCGGAGTCTGTC");
    link57.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 11, false, "TCATTCCACCGGAGTCTGTC");
    link57.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 36, false, "TCATTCCACCGGAGTCTGTC");
    if (linkset.find(link57) != linkset.end()) { ++linkset.at(link57); } else { linkset[link57] = 1; }
    auto link58 = Link();
    link58.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 79, false, "CAGACAGGCAGTCACATGAC");
    link58.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 81, false, "CAGACAGGCAGTCACATGAC");
    link58.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 79, false, "CAGACAGGCAGTCACATGAC");
    link58.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 81, false, "CAGACAGGCAGTCACATGAC");
    if (linkset.find(link58) != linkset.end()) { ++linkset.at(link58); } else { linkset[link58] = 1; }
    auto link59 = Link();
    link59.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 48, false, "AGTCTGTCTCATACCCAACC");
    link59.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 48, false, "AGTCTGTCTCATACCCAACC");
    if (linkset.find(link59) != linkset.end()) { ++linkset.at(link59); } else { linkset[link59] = 1; }
}

inline void fillLinksetTestdata_part3(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link60 = Link();
    link60.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 70, false, "CAGCAGGCACAGACAGGCAG");
    link60.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 72, false, "CAGCAGGCACAGACAGGCAG");
    link60.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 70, false, "CAGCAGGCACAGACAGGCAG");
    link60.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 72, false, "CAGCAGGCACAGACAGGCAG");
    if (linkset.find(link60) != linkset.end()) { ++linkset.at(link60); } else { linkset[link60] = 1; }
    auto link61 = Link();
    link61.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 52, false, "TGTCTCATACCCAACCAGAT");
    link61.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 52, false, "TGTCTCATACCCAACCAGAT");
    if (linkset.find(link61) != linkset.end()) { ++linkset.at(link61); } else { linkset[link61] = 1; }
    auto link62 = Link();
    link62.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 56, false, "TCATACCCAACCAGATTTCA");
    link62.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 56, false, "TCATACCCAACCAGATTTCA");
    if (linkset.find(link62) != linkset.end()) { ++linkset.at(link62); } else { linkset[link62] = 1; }
    auto link63 = Link();
    link63.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 65, false, "ACCAGATTTCAGTGGAGTGA");
    link63.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 39, false, "ACCAGATTTCAGTGGAGTGA");
    link63.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 65, false, "ACCAGATTTCAGTGGAGTGA");
    if (linkset.find(link63) != linkset.end()) { ++linkset.at(link63); } else { linkset[link63] = 1; }
    auto link64 = Link();
    link64.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 98, false, "AGTGGAGTGAAGTTCAGGAG");
    link64.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 97, false, "AGTGGAGTGAAGTTCAGGAG");
    if (linkset.find(link64) != linkset.end()) { ++linkset.at(link64); } else { linkset[link64] = 1; }
    auto link65 = Link();
    link65.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 17, false, "CATGTGCTTCTCTTGTCCTT");
    link65.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 17, false, "CATGTGCTTCTCTTGTCCTT");
    if (linkset.find(link65) != linkset.end()) { ++linkset.at(link65); } else { linkset[link65] = 1; }
    auto link66 = Link();
    link66.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 86, false, "GCAGTCACATGACAACCCAG");
    link66.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 88, false, "GCAGTCACATGACAACCCAG");
    link66.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 86, false, "GCAGTCACATGACAACCCAG");
    link66.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 88, false, "GCAGTCACATGACAACCCAG");
    if (linkset.find(link66) != linkset.end()) { ++linkset.at(link66); } else { linkset[link66] = 1; }
    auto link67 = Link();
    link67.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 23, false, "GTGTCTGCCTGTCTACACTT");
    link67.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 25, false, "GTGTCTGCCTGTCTACACTT");
    link67.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 23, false, "GTGTCTGCCTGTCTACACTT");
    link67.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 25, false, "GTGTCTGCCTGTCTACACTT");
    if (linkset.find(link67) != linkset.end()) { ++linkset.at(link67); } else { linkset[link67] = 1; }
    auto link68 = Link();
    link68.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 78, false, "ACAGACAGGCAGTCACATGA");
    link68.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 80, false, "ACAGACAGGCAGTCACATGA");
    link68.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 78, false, "ACAGACAGGCAGTCACATGA");
    link68.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 80, false, "ACAGACAGGCAGTCACATGA");
    if (linkset.find(link68) != linkset.end()) { ++linkset.at(link68); } else { linkset[link68] = 1; }
    auto link69 = Link();
    link69.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 54, false, "TCTCATACCCAACCAGATTT");
    link69.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 54, false, "TCTCATACCCAACCAGATTT");
    if (linkset.find(link69) != linkset.end()) { ++linkset.at(link69); } else { linkset[link69] = 1; }
    auto link70 = Link();
    link70.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 104, false, "GTGAAGTTCAGGAGGCATGG");
    link70.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 103, false, "GTGAAGTTCAGGAGGCATGG");
    if (linkset.find(link70) != linkset.end()) { ++linkset.at(link70); } else { linkset[link70] = 1; }
    auto link71 = Link();
    link71.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 76, false, "GCACAGACAGGCAGTCACAT");
    link71.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 78, false, "GCACAGACAGGCAGTCACAT");
    link71.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 76, false, "GCACAGACAGGCAGTCACAT");
    link71.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 78, false, "GCACAGACAGGCAGTCACAT");
    if (linkset.find(link71) != linkset.end()) { ++linkset.at(link71); } else { linkset[link71] = 1; }
    auto link72 = Link();
    link72.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 17, false, "TGTCATGTGTCTGCCTGTCT");
    link72.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 19, false, "TGTCATGTGTCTGCCTGTCT");
    link72.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 17, false, "TGTCATGTGTCTGCCTGTCT");
    link72.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 19, false, "TGTCATGTGTCTGCCTGTCT");
    if (linkset.find(link72) != linkset.end()) { ++linkset.at(link72); } else { linkset[link72] = 1; }
    auto link73 = Link();
    link73.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 35, false, "CTACACTTGCTGTGCAGAAC");
    link73.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 37, false, "CTACACTTGCTGTGCAGAAC");
    link73.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 35, false, "CTACACTTGCTGTGCAGAAC");
    link73.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 37, false, "CTACACTTGCTGTGCAGAAC");
    if (linkset.find(link73) != linkset.end()) { ++linkset.at(link73); } else { linkset[link73] = 1; }
    auto link74 = Link();
    link74.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 3, false, "TGGCTGGACAGAGTTGTCAT");
    link74.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 5, false, "TGGCTGGACAGAGTTGTCAT");
    link74.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 3, false, "TGGCTGGACAGAGTTGTCAT");
    link74.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 5, false, "TGGCTGGACAGAGTTGTCAT");
    if (linkset.find(link74) != linkset.end()) { ++linkset.at(link74); } else { linkset[link74] = 1; }
    auto link75 = Link();
    link75.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 63, false, "ACCTGTACAGCAGGCACAGA");
    link75.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 65, false, "ACCTGTACAGCAGGCACAGA");
    link75.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 63, false, "ACCTGTACAGCAGGCACAGA");
    link75.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 65, false, "ACCTGTACAGCAGGCACAGA");
    if (linkset.find(link75) != linkset.end()) { ++linkset.at(link75); } else { linkset[link75] = 1; }
    auto link76 = Link();
    link76.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 46, false, "GGAGTCTGTCTCATACCCAA");
    link76.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 46, false, "GGAGTCTGTCTCATACCCAA");
    if (linkset.find(link76) != linkset.end()) { ++linkset.at(link76); } else { linkset[link76] = 1; }
    auto link77 = Link();
    link77.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 40, false, "TCCACCGGAGTCTGTCTCAT");
    link77.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 40, false, "TCCACCGGAGTCTGTCTCAT");
    if (linkset.find(link77) != linkset.end()) { ++linkset.at(link77); } else { linkset[link77] = 1; }
    auto link78 = Link();
    link78.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 52, false, "AACATCCGCTCACCTGTACA");
    link78.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 54, false, "AACATCCGCTCACCTGTACA");
    link78.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 52, false, "AACATCCGCTCACCTGTACA");
    link78.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 54, false, "AACATCCGCTCACCTGTACA");
    if (linkset.find(link78) != linkset.end()) { ++linkset.at(link78); } else { linkset[link78] = 1; }
    auto link79 = Link();
    link79.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 1, false, "GGCCTGGCTGGACAGAGTTG");
    link79.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 1, false, "GGCCTGGCTGGACAGAGTTG");
    if (linkset.find(link79) != linkset.end()) { ++linkset.at(link79); } else { linkset[link79] = 1; }
}

inline void fillLinksetTestdata_part4(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link80 = Link();
    link80.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 45, false, "TGTGCAGAACATCCGCTCAC");
    link80.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 47, false, "TGTGCAGAACATCCGCTCAC");
    link80.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 45, false, "TGTGCAGAACATCCGCTCAC");
    link80.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 47, false, "TGTGCAGAACATCCGCTCAC");
    if (linkset.find(link80) != linkset.end()) { ++linkset.at(link80); } else { linkset[link80] = 1; }
    auto link81 = Link();
    link81.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 35, false, "TTCATTCCACCGGAGTCTGT");
    link81.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 10, false, "TTCATTCCACCGGAGTCTGT");
    link81.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 10, false, "TTCATTCCACCGGAGTCTGT");
    link81.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 35, false, "TTCATTCCACCGGAGTCTGT");
    if (linkset.find(link81) != linkset.end()) { ++linkset.at(link81); } else { linkset[link81] = 1; }
    auto link82 = Link();
    link82.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 65, false, "CTGTACAGCAGGCACAGACA");
    link82.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 67, false, "CTGTACAGCAGGCACAGACA");
    link82.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 65, false, "CTGTACAGCAGGCACAGACA");
    link82.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 67, false, "CTGTACAGCAGGCACAGACA");
    if (linkset.find(link82) != linkset.end()) { ++linkset.at(link82); } else { linkset[link82] = 1; }
    auto link83 = Link();
    link83.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 56, false, "TCCGCTCACCTGTACAGCAG");
    link83.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 58, false, "TCCGCTCACCTGTACAGCAG");
    link83.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 56, false, "TCCGCTCACCTGTACAGCAG");
    link83.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 58, false, "TCCGCTCACCTGTACAGCAG");
    if (linkset.find(link83) != linkset.end()) { ++linkset.at(link83); } else { linkset[link83] = 1; }
    auto link84 = Link();
    link84.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 37, false, "CATTCCACCGGAGTCTGTCT");
    link84.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 12, false, "CATTCCACCGGAGTCTGTCT");
    link84.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 12, false, "CATTCCACCGGAGTCTGTCT");
    link84.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 37, false, "CATTCCACCGGAGTCTGTCT");
    if (linkset.find(link84) != linkset.end()) { ++linkset.at(link84); } else { linkset[link84] = 1; }
    auto link85 = Link();
    link85.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 33, false, "GTCTACACTTGCTGTGCAGA");
    link85.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 35, false, "GTCTACACTTGCTGTGCAGA");
    link85.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 33, false, "GTCTACACTTGCTGTGCAGA");
    link85.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 35, false, "GTCTACACTTGCTGTGCAGA");
    if (linkset.find(link85) != linkset.end()) { ++linkset.at(link85); } else { linkset[link85] = 1; }
    auto link86 = Link();
    link86.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 34, false, "CTTCATTCCACCGGAGTCTG");
    link86.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 9, false, "CTTCATTCCACCGGAGTCTG");
    link86.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 9, false, "CTTCATTCCACCGGAGTCTG");
    link86.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 34, false, "CTTCATTCCACCGGAGTCTG");
    if (linkset.find(link86) != linkset.end()) { ++linkset.at(link86); } else { linkset[link86] = 1; }
    auto link87 = Link();
    link87.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 62, false, "CCAACCAGATTTCAGTGGAG");
    link87.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 36, false, "CCAACCAGATTTCAGTGGAG");
    link87.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 62, false, "CCAACCAGATTTCAGTGGAG");
    if (linkset.find(link87) != linkset.end()) { ++linkset.at(link87); } else { linkset[link87] = 1; }
    auto link88 = Link();
    link88.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 16, false, "TTGTCATGTGTCTGCCTGTC");
    link88.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 18, false, "TTGTCATGTGTCTGCCTGTC");
    link88.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 16, false, "TTGTCATGTGTCTGCCTGTC");
    link88.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 18, false, "TTGTCATGTGTCTGCCTGTC");
    if (linkset.find(link88) != linkset.end()) { ++linkset.at(link88); } else { linkset[link88] = 1; }
    auto link89 = Link();
    link89.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 53, false, "GTCTCATACCCAACCAGATT");
    link89.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 53, false, "GTCTCATACCCAACCAGATT");
    if (linkset.find(link89) != linkset.end()) { ++linkset.at(link89); } else { linkset[link89] = 1; }
    auto link90 = Link();
    link90.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 99, false, "GTGGAGTGAAGTTCAGGAGG");
    link90.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 98, false, "GTGGAGTGAAGTTCAGGAGG");
    if (linkset.find(link90) != linkset.end()) { ++linkset.at(link90); } else { linkset[link90] = 1; }
    auto link91 = Link();
    link91.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 55, false, "CTCATACCCAACCAGATTTC");
    link91.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 55, false, "CTCATACCCAACCAGATTTC");
    if (linkset.find(link91) != linkset.end()) { ++linkset.at(link91); } else { linkset[link91] = 1; }
    auto link92 = Link();
    link92.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 106, false, "GAAGTTCAGGAGGCATGGAG");
    link92.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 105, false, "GAAGTTCAGGAGGCATGGAG");
    if (linkset.find(link92) != linkset.end()) { ++linkset.at(link92); } else { linkset[link92] = 1; }
    auto link93 = Link();
    link93.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 48, false, "GCAGAACATCCGCTCACCTG");
    link93.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 50, false, "GCAGAACATCCGCTCACCTG");
    link93.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 48, false, "GCAGAACATCCGCTCACCTG");
    link93.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 50, false, "GCAGAACATCCGCTCACCTG");
    if (linkset.find(link93) != linkset.end()) { ++linkset.at(link93); } else { linkset[link93] = 1; }
    auto link94 = Link();
    link94.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 55, false, "ATCCGCTCACCTGTACAGCA");
    link94.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 57, false, "ATCCGCTCACCTGTACAGCA");
    link94.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 55, false, "ATCCGCTCACCTGTACAGCA");
    link94.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 57, false, "ATCCGCTCACCTGTACAGCA");
    if (linkset.find(link94) != linkset.end()) { ++linkset.at(link94); } else { linkset[link94] = 1; }
    auto link95 = Link();
    link95.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 101, false, "GGAGTGAAGTTCAGGAGGCA");
    link95.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 100, false, "GGAGTGAAGTTCAGGAGGCA");
    if (linkset.find(link95) != linkset.end()) { ++linkset.at(link95); } else { linkset[link95] = 1; }
    auto link96 = Link();
    link96.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 37, false, "ACACTTGCTGTGCAGAACAT");
    link96.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 39, false, "ACACTTGCTGTGCAGAACAT");
    link96.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 37, false, "ACACTTGCTGTGCAGAACAT");
    link96.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 39, false, "ACACTTGCTGTGCAGAACAT");
    if (linkset.find(link96) != linkset.end()) { ++linkset.at(link96); } else { linkset[link96] = 1; }
    auto link97 = Link();
    link97.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 108, false, "AGTTCAGGAGGCATGGAGCT");
    link97.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 107, false, "AGTTCAGGAGGCATGGAGCT");
    if (linkset.find(link97) != linkset.end()) { ++linkset.at(link97); } else { linkset[link97] = 1; }
    auto link98 = Link();
    link98.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 39, false, "ACTTGCTGTGCAGAACATCC");
    link98.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 41, false, "ACTTGCTGTGCAGAACATCC");
    link98.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 39, false, "ACTTGCTGTGCAGAACATCC");
    link98.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 41, false, "ACTTGCTGTGCAGAACATCC");
    if (linkset.find(link98) != linkset.end()) { ++linkset.at(link98); } else { linkset[link98] = 1; }
}

inline void fillLinksetTestdata_part5(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link99 = Link();
    link99.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 84, false, "AGGCAGTCACATGACAACCC");
    link99.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 86, false, "AGGCAGTCACATGACAACCC");
    link99.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 84, false, "AGGCAGTCACATGACAACCC");
    link99.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 86, false, "AGGCAGTCACATGACAACCC");
    if (linkset.find(link99) != linkset.end()) { ++linkset.at(link99); } else { linkset[link99] = 1; }
    auto link100 = Link();
    link100.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 47, false, "GAGTCTGTCTCATACCCAAC");
    link100.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 47, false, "GAGTCTGTCTCATACCCAAC");
    if (linkset.find(link100) != linkset.end()) { ++linkset.at(link100); } else { linkset[link100] = 1; }
    auto link101 = Link();
    link101.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 63, false, "CAACCAGATTTCAGTGGAGT");
    link101.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 37, false, "CAACCAGATTTCAGTGGAGT");
    link101.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 63, false, "CAACCAGATTTCAGTGGAGT");
    if (linkset.find(link101) != linkset.end()) { ++linkset.at(link101); } else { linkset[link101] = 1; }
    auto link102 = Link();
    link102.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 75, false, "GGCACAGACAGGCAGTCACA");
    link102.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 77, false, "GGCACAGACAGGCAGTCACA");
    link102.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 75, false, "GGCACAGACAGGCAGTCACA");
    link102.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 77, false, "GGCACAGACAGGCAGTCACA");
    if (linkset.find(link102) != linkset.end()) { ++linkset.at(link102); } else { linkset[link102] = 1; }
    auto link103 = Link();
    link103.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 81, false, "GACAGGCAGTCACATGACAA");
    link103.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 83, false, "GACAGGCAGTCACATGACAA");
    link103.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 81, false, "GACAGGCAGTCACATGACAA");
    link103.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 83, false, "GACAGGCAGTCACATGACAA");
    if (linkset.find(link103) != linkset.end()) { ++linkset.at(link103); } else { linkset[link103] = 1; }
    auto link104 = Link();
    link104.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 2, false, "CTGGCTGGACAGAGTTGTCA");
    link104.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 4, false, "CTGGCTGGACAGAGTTGTCA");
    link104.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 2, false, "CTGGCTGGACAGAGTTGTCA");
    link104.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 4, false, "CTGGCTGGACAGAGTTGTCA");
    if (linkset.find(link104) != linkset.end()) { ++linkset.at(link104); } else { linkset[link104] = 1; }
    auto link105 = Link();
    link105.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 33, false, "CCTTCATTCCACCGGAGTCT");
    link105.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 8, false, "CCTTCATTCCACCGGAGTCT");
    link105.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 8, false, "CCTTCATTCCACCGGAGTCT");
    link105.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 33, false, "CCTTCATTCCACCGGAGTCT");
    if (linkset.find(link105) != linkset.end()) { ++linkset.at(link105); } else { linkset[link105] = 1; }
    auto link106 = Link();
    link106.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 58, false, "CGCTCACCTGTACAGCAGGC");
    link106.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 60, false, "CGCTCACCTGTACAGCAGGC");
    link106.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 58, false, "CGCTCACCTGTACAGCAGGC");
    link106.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 60, false, "CGCTCACCTGTACAGCAGGC");
    if (linkset.find(link106) != linkset.end()) { ++linkset.at(link106); } else { linkset[link106] = 1; }
    auto link107 = Link();
    link107.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 6, false, "CTGGACAGAGTTGTCATGTG");
    link107.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 8, false, "CTGGACAGAGTTGTCATGTG");
    link107.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 6, false, "CTGGACAGAGTTGTCATGTG");
    link107.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 8, false, "CTGGACAGAGTTGTCATGTG");
    if (linkset.find(link107) != linkset.end()) { ++linkset.at(link107); } else { linkset[link107] = 1; }
    auto link108 = Link();
    link108.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 24, false, "TTCTCTTGTCCTTCATTCCA");
    link108.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 24, false, "TTCTCTTGTCCTTCATTCCA");
    if (linkset.find(link108) != linkset.end()) { ++linkset.at(link108); } else { linkset[link108] = 1; }
    auto link109 = Link();
    link109.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 42, false, "TGCTGTGCAGAACATCCGCT");
    link109.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 44, false, "TGCTGTGCAGAACATCCGCT");
    link109.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 42, false, "TGCTGTGCAGAACATCCGCT");
    link109.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 44, false, "TGCTGTGCAGAACATCCGCT");
    if (linkset.find(link109) != linkset.end()) { ++linkset.at(link109); } else { linkset[link109] = 1; }
    auto link110 = Link();
    link110.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 43, false, "ACCGGAGTCTGTCTCATACC");
    link110.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 43, false, "ACCGGAGTCTGTCTCATACC");
    if (linkset.find(link110) != linkset.end()) { ++linkset.at(link110); } else { linkset[link110] = 1; }
    auto link111 = Link();
    link111.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 21, false, "ATGTGTCTGCCTGTCTACAC");
    link111.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 23, false, "ATGTGTCTGCCTGTCTACAC");
    link111.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 21, false, "ATGTGTCTGCCTGTCTACAC");
    link111.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 23, false, "ATGTGTCTGCCTGTCTACAC");
    if (linkset.find(link111) != linkset.end()) { ++linkset.at(link111); } else { linkset[link111] = 1; }
    auto link112 = Link();
    link112.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 68, false, "TACAGCAGGCACAGACAGGC");
    link112.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 70, false, "TACAGCAGGCACAGACAGGC");
    link112.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 68, false, "TACAGCAGGCACAGACAGGC");
    link112.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 70, false, "TACAGCAGGCACAGACAGGC");
    if (linkset.find(link112) != linkset.end()) { ++linkset.at(link112); } else { linkset[link112] = 1; }
    auto link113 = Link();
    link113.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 62, false, "CACCTGTACAGCAGGCACAG");
    link113.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 64, false, "CACCTGTACAGCAGGCACAG");
    link113.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 62, false, "CACCTGTACAGCAGGCACAG");
    link113.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 64, false, "CACCTGTACAGCAGGCACAG");
    if (linkset.find(link113) != linkset.end()) { ++linkset.at(link113); } else { linkset[link113] = 1; }
    auto link114 = Link();
    link114.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 73, false, "CAGGCACAGACAGGCAGTCA");
    link114.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 75, false, "CAGGCACAGACAGGCAGTCA");
    link114.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 73, false, "CAGGCACAGACAGGCAGTCA");
    link114.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 75, false, "CAGGCACAGACAGGCAGTCA");
    if (linkset.find(link114) != linkset.end()) { ++linkset.at(link114); } else { linkset[link114] = 1; }
    auto link115 = Link();
    link115.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 94, false, "TTTCAGTGGAGTGAAGTTCA");
    link115.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 93, false, "TTTCAGTGGAGTGAAGTTCA");
    if (linkset.find(link115) != linkset.end()) { ++linkset.at(link115); } else { linkset[link115] = 1; }
    auto link116 = Link();
    link116.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 85, false, "GGCAGTCACATGACAACCCA");
    link116.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 87, false, "GGCAGTCACATGACAACCCA");
    link116.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 85, false, "GGCAGTCACATGACAACCCA");
    link116.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 87, false, "GGCAGTCACATGACAACCCA");
    if (linkset.find(link116) != linkset.end()) { ++linkset.at(link116); } else { linkset[link116] = 1; }
    auto link117 = Link();
    link117.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 49, false, "CAGAACATCCGCTCACCTGT");
    link117.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 51, false, "CAGAACATCCGCTCACCTGT");
    link117.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 49, false, "CAGAACATCCGCTCACCTGT");
    link117.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 51, false, "CAGAACATCCGCTCACCTGT");
    if (linkset.find(link117) != linkset.end()) { ++linkset.at(link117); } else { linkset[link117] = 1; }
}

inline void fillLinksetTestdata_part6(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link118 = Link();
    link118.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 14, false, "ATCCATGTGCTTCTCTTGTC");
    link118.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 14, false, "ATCCATGTGCTTCTCTTGTC");
    if (linkset.find(link118) != linkset.end()) { ++linkset.at(link118); } else { linkset[link118] = 1; }
    auto link119 = Link();
    link119.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 31, false, "GTCCTTCATTCCACCGGAGT");
    link119.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 6, false, "GTCCTTCATTCCACCGGAGT");
    link119.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 6, false, "GTCCTTCATTCCACCGGAGT");
    link119.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 31, false, "GTCCTTCATTCCACCGGAGT");
    if (linkset.find(link119) != linkset.end()) { ++linkset.at(link119); } else { linkset[link119] = 1; }
    auto link120 = Link();
    link120.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 61, false, "TCACCTGTACAGCAGGCACA");
    link120.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 63, false, "TCACCTGTACAGCAGGCACA");
    link120.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 61, false, "TCACCTGTACAGCAGGCACA");
    link120.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 63, false, "TCACCTGTACAGCAGGCACA");
    if (linkset.find(link120) != linkset.end()) { ++linkset.at(link120); } else { linkset[link120] = 1; }
    auto link121 = Link();
    link121.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 19, false, "TCATGTGTCTGCCTGTCTAC");
    link121.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 21, false, "TCATGTGTCTGCCTGTCTAC");
    link121.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 19, false, "TCATGTGTCTGCCTGTCTAC");
    link121.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 21, false, "TCATGTGTCTGCCTGTCTAC");
    if (linkset.find(link121) != linkset.end()) { ++linkset.at(link121); } else { linkset[link121] = 1; }
    auto link122 = Link();
    link122.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 29, false, "GCCTGTCTACACTTGCTGTG");
    link122.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 31, false, "GCCTGTCTACACTTGCTGTG");
    link122.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 29, false, "GCCTGTCTACACTTGCTGTG");
    link122.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 31, false, "GCCTGTCTACACTTGCTGTG");
    if (linkset.find(link122) != linkset.end()) { ++linkset.at(link122); } else { linkset[link122] = 1; }
    auto link123 = Link();
    link123.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 67, false, "CAGATTTCAGTGGAGTGAAG");
    link123.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 41, false, "CAGATTTCAGTGGAGTGAAG");
    link123.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 67, false, "CAGATTTCAGTGGAGTGAAG");
    if (linkset.find(link123) != linkset.end()) { ++linkset.at(link123); } else { linkset[link123] = 1; }
    auto link124 = Link();
    link124.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 38, false, "ATTCCACCGGAGTCTGTCTC");
    link124.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 38, false, "ATTCCACCGGAGTCTGTCTC");
    if (linkset.find(link124) != linkset.end()) { ++linkset.at(link124); } else { linkset[link124] = 1; }
    auto link125 = Link();
    link125.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 69, false, "ACAGCAGGCACAGACAGGCA");
    link125.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 71, false, "ACAGCAGGCACAGACAGGCA");
    link125.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 69, false, "ACAGCAGGCACAGACAGGCA");
    link125.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 71, false, "ACAGCAGGCACAGACAGGCA");
    if (linkset.find(link125) != linkset.end()) { ++linkset.at(link125); } else { linkset[link125] = 1; }
    auto link126 = Link();
    link126.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 57, false, "CCGCTCACCTGTACAGCAGG");
    link126.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 59, false, "CCGCTCACCTGTACAGCAGG");
    link126.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 57, false, "CCGCTCACCTGTACAGCAGG");
    link126.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 59, false, "CCGCTCACCTGTACAGCAGG");
    if (linkset.find(link126) != linkset.end()) { ++linkset.at(link126); } else { linkset[link126] = 1; }
    auto link127 = Link();
    link127.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 64, false, "CCTGTACAGCAGGCACAGAC");
    link127.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 66, false, "CCTGTACAGCAGGCACAGAC");
    link127.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 64, false, "CCTGTACAGCAGGCACAGAC");
    link127.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 66, false, "CCTGTACAGCAGGCACAGAC");
    if (linkset.find(link127) != linkset.end()) { ++linkset.at(link127); } else { linkset[link127] = 1; }
    auto link128 = Link();
    link128.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 43, false, "GCTGTGCAGAACATCCGCTC");
    link128.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 45, false, "GCTGTGCAGAACATCCGCTC");
    link128.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 43, false, "GCTGTGCAGAACATCCGCTC");
    link128.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 45, false, "GCTGTGCAGAACATCCGCTC");
    if (linkset.find(link128) != linkset.end()) { ++linkset.at(link128); } else { linkset[link128] = 1; }
    auto link129 = Link();
    link129.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 41, false, "TTGCTGTGCAGAACATCCGC");
    link129.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 43, false, "TTGCTGTGCAGAACATCCGC");
    link129.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 41, false, "TTGCTGTGCAGAACATCCGC");
    link129.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 43, false, "TTGCTGTGCAGAACATCCGC");
    if (linkset.find(link129) != linkset.end()) { ++linkset.at(link129); } else { linkset[link129] = 1; }
    auto link130 = Link();
    link130.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 44, false, "CCGGAGTCTGTCTCATACCC");
    link130.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 44, false, "CCGGAGTCTGTCTCATACCC");
    if (linkset.find(link130) != linkset.end()) { ++linkset.at(link130); } else { linkset[link130] = 1; }
    auto link131 = Link();
    link131.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 42, false, "CACCGGAGTCTGTCTCATAC");
    link131.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 42, false, "CACCGGAGTCTGTCTCATAC");
    if (linkset.find(link131) != linkset.end()) { ++linkset.at(link131); } else { linkset[link131] = 1; }
    auto link132 = Link();
    link132.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 7, false, "TGGACAGAGTTGTCATGTGT");
    link132.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 9, false, "TGGACAGAGTTGTCATGTGT");
    link132.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 7, false, "TGGACAGAGTTGTCATGTGT");
    link132.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 9, false, "TGGACAGAGTTGTCATGTGT");
    if (linkset.find(link132) != linkset.end()) { ++linkset.at(link132); } else { linkset[link132] = 1; }
    auto link133 = Link();
    link133.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 60, false, "ACCCAACCAGATTTCAGTGG");
    link133.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 60, false, "ACCCAACCAGATTTCAGTGG");
    if (linkset.find(link133) != linkset.end()) { ++linkset.at(link133); } else { linkset[link133] = 1; }
    auto link134 = Link();
    link134.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 26, false, "TCTGCCTGTCTACACTTGCT");
    link134.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 28, false, "TCTGCCTGTCTACACTTGCT");
    link134.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 26, false, "TCTGCCTGTCTACACTTGCT");
    link134.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 28, false, "TCTGCCTGTCTACACTTGCT");
    if (linkset.find(link134) != linkset.end()) { ++linkset.at(link134); } else { linkset[link134] = 1; }
    auto link135 = Link();
    link135.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 12, false, "CAATCCATGTGCTTCTCTTG");
    link135.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 12, false, "CAATCCATGTGCTTCTCTTG");
    if (linkset.find(link135) != linkset.end()) { ++linkset.at(link135); } else { linkset[link135] = 1; }
    auto link136 = Link();
    link136.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 66, false, "CCAGATTTCAGTGGAGTGAA");
    link136.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 40, false, "CCAGATTTCAGTGGAGTGAA");
    link136.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 66, false, "CCAGATTTCAGTGGAGTGAA");
    if (linkset.find(link136) != linkset.end()) { ++linkset.at(link136); } else { linkset[link136] = 1; }
}

inline void fillLinksetTestdata_part7(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link137 = Link();
    link137.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 107, false, "AAGTTCAGGAGGCATGGAGC");
    link137.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 106, false, "AAGTTCAGGAGGCATGGAGC");
    if (linkset.find(link137) != linkset.end()) { ++linkset.at(link137); } else { linkset[link137] = 1; }
    auto link138 = Link();
    link138.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 59, false, "GCTCACCTGTACAGCAGGCA");
    link138.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 61, false, "GCTCACCTGTACAGCAGGCA");
    link138.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 59, false, "GCTCACCTGTACAGCAGGCA");
    link138.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 61, false, "GCTCACCTGTACAGCAGGCA");
    if (linkset.find(link138) != linkset.end()) { ++linkset.at(link138); } else { linkset[link138] = 1; }
    auto link139 = Link();
    link139.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 28, false, "TGCCTGTCTACACTTGCTGT");
    link139.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 30, false, "TGCCTGTCTACACTTGCTGT");
    link139.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 28, false, "TGCCTGTCTACACTTGCTGT");
    link139.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 30, false, "TGCCTGTCTACACTTGCTGT");
    if (linkset.find(link139) != linkset.end()) { ++linkset.at(link139); } else { linkset[link139] = 1; }
    auto link140 = Link();
    link140.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 80, false, "AGACAGGCAGTCACATGACA");
    link140.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 82, false, "AGACAGGCAGTCACATGACA");
    link140.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 80, false, "AGACAGGCAGTCACATGACA");
    link140.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 82, false, "AGACAGGCAGTCACATGACA");
    if (linkset.find(link140) != linkset.end()) { ++linkset.at(link140); } else { linkset[link140] = 1; }
    auto link141 = Link();
    link141.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 44, false, "CTGTGCAGAACATCCGCTCA");
    link141.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 46, false, "CTGTGCAGAACATCCGCTCA");
    link141.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 44, false, "CTGTGCAGAACATCCGCTCA");
    link141.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 46, false, "CTGTGCAGAACATCCGCTCA");
    if (linkset.find(link141) != linkset.end()) { ++linkset.at(link141); } else { linkset[link141] = 1; }
    auto link142 = Link();
    link142.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 46, false, "GTGCAGAACATCCGCTCACC");
    link142.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 48, false, "GTGCAGAACATCCGCTCACC");
    link142.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 46, false, "GTGCAGAACATCCGCTCACC");
    link142.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 48, false, "GTGCAGAACATCCGCTCACC");
    if (linkset.find(link142) != linkset.end()) { ++linkset.at(link142); } else { linkset[link142] = 1; }
    auto link143 = Link();
    link143.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 34, false, "TCTACACTTGCTGTGCAGAA");
    link143.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 36, false, "TCTACACTTGCTGTGCAGAA");
    link143.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 34, false, "TCTACACTTGCTGTGCAGAA");
    link143.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 36, false, "TCTACACTTGCTGTGCAGAA");
    if (linkset.find(link143) != linkset.end()) { ++linkset.at(link143); } else { linkset[link143] = 1; }
    auto link144 = Link();
    link144.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 22, false, "TGTGTCTGCCTGTCTACACT");
    link144.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 24, false, "TGTGTCTGCCTGTCTACACT");
    link144.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 22, false, "TGTGTCTGCCTGTCTACACT");
    link144.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 24, false, "TGTGTCTGCCTGTCTACACT");
    if (linkset.find(link144) != linkset.end()) { ++linkset.at(link144); } else { linkset[link144] = 1; }
    auto link145 = Link();
    link145.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 15, false, "TCCATGTGCTTCTCTTGTCC");
    link145.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 15, false, "TCCATGTGCTTCTCTTGTCC");
    if (linkset.find(link145) != linkset.end()) { ++linkset.at(link145); } else { linkset[link145] = 1; }
    auto link146 = Link();
    link146.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 14, false, "AGTTGTCATGTGTCTGCCTG");
    link146.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 16, false, "AGTTGTCATGTGTCTGCCTG");
    link146.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 14, false, "AGTTGTCATGTGTCTGCCTG");
    link146.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 16, false, "AGTTGTCATGTGTCTGCCTG");
    if (linkset.find(link146) != linkset.end()) { ++linkset.at(link146); } else { linkset[link146] = 1; }
    auto link147 = Link();
    link147.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 66, false, "TGTACAGCAGGCACAGACAG");
    link147.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 68, false, "TGTACAGCAGGCACAGACAG");
    link147.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 66, false, "TGTACAGCAGGCACAGACAG");
    link147.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 68, false, "TGTACAGCAGGCACAGACAG");
    if (linkset.find(link147) != linkset.end()) { ++linkset.at(link147); } else { linkset[link147] = 1; }
    auto link148 = Link();
    link148.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 67, false, "GTACAGCAGGCACAGACAGG");
    link148.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 69, false, "GTACAGCAGGCACAGACAGG");
    link148.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 67, false, "GTACAGCAGGCACAGACAGG");
    link148.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 69, false, "GTACAGCAGGCACAGACAGG");
    if (linkset.find(link148) != linkset.end()) { ++linkset.at(link148); } else { linkset[link148] = 1; }
    auto link149 = Link();
    link149.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 4, false, "GGCTGGACAGAGTTGTCATG");
    link149.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 6, false, "GGCTGGACAGAGTTGTCATG");
    link149.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 4, false, "GGCTGGACAGAGTTGTCATG");
    link149.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 6, false, "GGCTGGACAGAGTTGTCATG");
    if (linkset.find(link149) != linkset.end()) { ++linkset.at(link149); } else { linkset[link149] = 1; }
    auto link150 = Link();
    link150.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 30, false, "TGTCCTTCATTCCACCGGAG");
    link150.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 5, false, "TGTCCTTCATTCCACCGGAG");
    link150.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 5, false, "TGTCCTTCATTCCACCGGAG");
    link150.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 30, false, "TGTCCTTCATTCCACCGGAG");
    if (linkset.find(link150) != linkset.end()) { ++linkset.at(link150); } else { linkset[link150] = 1; }
    auto link151 = Link();
    link151.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 88, false, "AGTCACATGACAACCCAGCC");
    link151.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 90, false, "AGTCACATGACAACCCAGCC");
    link151.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 88, false, "AGTCACATGACAACCCAGCC");
    link151.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 90, false, "AGTCACATGACAACCCAGCC");
    if (linkset.find(link151) != linkset.end()) { ++linkset.at(link151); } else { linkset[link151] = 1; }
    auto link152 = Link();
    link152.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 97, false, "CAGTGGAGTGAAGTTCAGGA");
    link152.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 96, false, "CAGTGGAGTGAAGTTCAGGA");
    if (linkset.find(link152) != linkset.end()) { ++linkset.at(link152); } else { linkset[link152] = 1; }
    auto link153 = Link();
    link153.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 13, false, "AATCCATGTGCTTCTCTTGT");
    link153.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 13, false, "AATCCATGTGCTTCTCTTGT");
    if (linkset.find(link153) != linkset.end()) { ++linkset.at(link153); } else { linkset[link153] = 1; }
    auto link154 = Link();
    link154.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 51, false, "GAACATCCGCTCACCTGTAC");
    link154.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 53, false, "GAACATCCGCTCACCTGTAC");
    link154.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 51, false, "GAACATCCGCTCACCTGTAC");
    link154.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 53, false, "GAACATCCGCTCACCTGTAC");
    if (linkset.find(link154) != linkset.end()) { ++linkset.at(link154); } else { linkset[link154] = 1; }
}

inline void fillLinksetTestdata_part8(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    auto link155 = Link();
    link155.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 74, false, "AGGCACAGACAGGCAGTCAC");
    link155.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 76, false, "AGGCACAGACAGGCAGTCAC");
    link155.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 74, false, "AGGCACAGACAGGCAGTCAC");
    link155.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 76, false, "AGGCACAGACAGGCAGTCAC");
    if (linkset.find(link155) != linkset.end()) { ++linkset.at(link155); } else { linkset[link155] = 1; }
    auto link156 = Link();
    link156.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 36, false, "TACACTTGCTGTGCAGAACA");
    link156.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 38, false, "TACACTTGCTGTGCAGAACA");
    link156.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 36, false, "TACACTTGCTGTGCAGAACA");
    link156.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 38, false, "TACACTTGCTGTGCAGAACA");
    if (linkset.find(link156) != linkset.end()) { ++linkset.at(link156); } else { linkset[link156] = 1; }
    auto link157 = Link();
    link157.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 50, false, "TCTGTCTCATACCCAACCAG");
    link157.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 50, false, "TCTGTCTCATACCCAACCAG");
    if (linkset.find(link157) != linkset.end()) { ++linkset.at(link157); } else { linkset[link157] = 1; }
    auto link158 = Link();
    link158.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 32, false, "TCCTTCATTCCACCGGAGTC");
    link158.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 7, false, "TCCTTCATTCCACCGGAGTC");
    link158.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|193507463|68|-|195471971|68")), 7, false, "TCCTTCATTCCACCGGAGTC");
    link158.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 32, false, "TCCTTCATTCCACCGGAGTC");
    if (linkset.find(link158) != linkset.end()) { ++linkset.at(link158); } else { linkset[link158] = 1; }
    auto link159 = Link();
    link159.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 11, false, "CAGAGTTGTCATGTGTCTGC");
    link159.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 13, false, "CAGAGTTGTCATGTGTCTGC");
    link159.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 11, false, "CAGAGTTGTCATGTGTCTGC");
    link159.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 13, false, "CAGAGTTGTCATGTGTCTGC");
    if (linkset.find(link159) != linkset.end()) { ++linkset.at(link159); } else { linkset[link159] = 1; }
    auto link160 = Link();
    link160.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 102, false, "GAGTGAAGTTCAGGAGGCAT");
    link160.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 101, false, "GAGTGAAGTTCAGGAGGCAT");
    if (linkset.find(link160) != linkset.end()) { ++linkset.at(link160); } else { linkset[link160] = 1; }
    auto link161 = Link();
    link161.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 25, false, "TCTCTTGTCCTTCATTCCAC");
    link161.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602120|2525816|67|-|10179728|67")), 0, false, "TCTCTTGTCCTTCATTCCAC");
    link161.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 25, false, "TCTCTTGTCCTTCATTCCAC");
    if (linkset.find(link161) != linkset.end()) { ++linkset.at(link161); } else { linkset[link161] = 1; }
    auto link162 = Link();
    link162.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 24, false, "TGTCTGCCTGTCTACACTTG");
    link162.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 26, false, "TGTCTGCCTGTCTACACTTG");
    link162.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 24, false, "TGTCTGCCTGTCTACACTTG");
    link162.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 26, false, "TGTCTGCCTGTCTACACTTG");
    if (linkset.find(link162) != linkset.end()) { ++linkset.at(link162); } else { linkset[link162] = 1; }
    auto link163 = Link();
    link163.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 96, false, "TCAGTGGAGTGAAGTTCAGG");
    link163.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 95, false, "TCAGTGGAGTGAAGTTCAGG");
    if (linkset.find(link163) != linkset.end()) { ++linkset.at(link163); } else { linkset[link163] = 1; }
    auto link164 = Link();
    link164.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("mm10_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("mm10|chr1|162223368|110|+|195471971|110")), 9, false, "GACAGAGTTGTCATGTGTCT");
    link164.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|172138798|110|-|248956422|110")), 11, false, "GACAGAGTTGTCATGTGTCT");
    link164.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hetGla2_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hetGla2|JH602084|9538844|110|+|20331017|110")), 9, false, "GACAGAGTTGTCATGTGTCT");
    link164.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|30111737|110|-|227556264|110")), 11, false, "GACAGAGTTGTCATGTGTCT");
    if (linkset.find(link164) != linkset.end()) { ++linkset.at(link164); } else { linkset[link164] = 1; }
    auto link165 = Link();
    link165.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("hg38_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("hg38|chr1|209432133|110|+|248956422|110")), 57, false, "CATACCCAACCAGATTTCAG");
    link165.insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("macFas5_orthologs")), static_cast<uint32_t>(idMap.querySequenceID_LEGACY("macFas5|chr1|67985998|109|+|227556264|109")), 57, false, "CATACCCAACCAGATTTCAG");
    if (linkset.find(link165) != linkset.end()) { ++linkset.at(link165); } else { linkset[link165] = 1; }
}

inline void fillLinksetTestdata(typename Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>::LinksetType & linkset, IdentifierMapping const & idMap) {
    fillLinksetTestdata_part0(linkset, idMap);
    fillLinksetTestdata_part1(linkset, idMap);
    fillLinksetTestdata_part2(linkset, idMap);
    fillLinksetTestdata_part3(linkset, idMap);
    fillLinksetTestdata_part4(linkset, idMap);
    fillLinksetTestdata_part5(linkset, idMap);
    fillLinksetTestdata_part6(linkset, idMap);
    fillLinksetTestdata_part7(linkset, idMap);
    fillLinksetTestdata_part8(linkset, idMap);
}

#endif // FILLLINKSETTESTDATA_H
