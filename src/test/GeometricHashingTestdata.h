#ifndef MANUALTESTSET_H
#define MANUALTESTSET_H

#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "Cube.h"
#include "Cubeset.h"
#include "IdentifierMapping.h"
#include "Link.h"

// "Variables declared constexpr or const, and whose value is fixed for the duration of the program, are named with a leading "k" followed by mixed case"
// https://google.github.io/styleguide/cppguide.html#Constant_Names
auto const kTestGraph = "./data/test/graph";
std::vector<std::string> const kTestSpecies{"species1",
                                            "species2",
                                            "species3",
                                            "species4",
                                            "species5"};
std::vector<std::string> const kTestSequences{"species1_sequence1",
                                              "species2_sequence1",
                                              "species2_sequence2",
                                              "species3_sequence1",
                                              "species3_sequence2",
                                              "species4_sequence1",
                                              "species4_sequence2",
                                              "species4_sequence3",
                                              "species5_sequence1",
                                              "species5_sequence2"};
double const kAmplitudeDelta = 2.; // fix
double const kAmplitudeLink = 2.;  // fix
double const kMidDelta = 0.;   // fix
double const kMidLink = 0.;    // fix
double const kSteepnessDelta = 0.5; // default?
double const kSteepnessLink = 1.;   // default?
double const kVerticalShiftDelta = -1.;    // fix
double const kVerticalShiftLink = -1.;     // fix

auto const kNumGenomes = static_cast<uint16_t>(kTestSpecies.size());

class ManualSets {
public:
    /* Manually create the Linkset, Cubeset, Hasse diagram
     * that are expected from the testdata
     * ~~~
     * compare this with data/test/createTestdata.R
     * ~~~
     * coord-binsize = 1, i.e. work with exact positions -> k=5-mers are
     *     concatenated by "N", so positions are 0, 6, 12, ...
     * ~~~
     * See Labbook 2018-11-13 */
    ManualSets(IdentifierMapping & idMap)
        : expectedCubes_{}, expectedLinkset_{}, expectedLinksInCubes_{},
          expectedPredecessors_{}, expectedPredecessorsComplete_{},
          expectedScore_{} {
        // get internal IDs
        auto spec0 = idMap.queryGenomeID(kTestSpecies[0]);
        auto spec1 = idMap.queryGenomeID(kTestSpecies[1]);
        auto spec2 = idMap.queryGenomeID(kTestSpecies[2]);
        auto spec3 = idMap.queryGenomeID(kTestSpecies[3]);
        auto spec4 = idMap.queryGenomeID(kTestSpecies[4]);
        auto seq00 = idMap.querySequenceID(kTestSequences[0]);
        auto seq10 = idMap.querySequenceID(kTestSequences[1]);
        auto seq11 = idMap.querySequenceID(kTestSequences[2]);
        auto seq20 = idMap.querySequenceID(kTestSequences[3]);
        auto seq21 = idMap.querySequenceID(kTestSequences[4]);
        auto seq30 = idMap.querySequenceID(kTestSequences[5]);
        auto seq31 = idMap.querySequenceID(kTestSequences[6]);
        auto seq32 = idMap.querySequenceID(kTestSequences[7]);
        auto seq40 = idMap.querySequenceID(kTestSequences[8]);
        auto seq41 = idMap.querySequenceID(kTestSequences[9]);

        // create Linkset
        auto link = std::make_shared<Link>();       // 00 - 10
        link->insertOccurrence(spec0, seq00, 0, false, "TTGGC");
        link->insertOccurrence(spec1, seq10, 0, false, "TTGGC");
        auto c00_10 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_10);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_10].insert(link);

        link = std::make_shared<Link>();       // 00 - 10
        link->insertOccurrence(spec0, seq00, 6, false, "TATGC");
        link->insertOccurrence(spec1, seq10, 6, false, "TATGC");
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_10].insert(link);
        expectedScore_[c00_10] = 4;


        link = std::make_shared<Link>();            // 00 - 11
        link->insertOccurrence(spec0, seq00, 12, false, "GTAAC");
        link->insertOccurrence(spec1, seq11, 0, false, "GTAAC");
        auto c00_11 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_11);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_11].insert(link);
        expectedScore_[c00_11] = 1;

        link = std::make_shared<Link>();            // 00 - 21
        link->insertOccurrence(spec0, seq00, 18, false, "TCAGG");
        link->insertOccurrence(spec2, seq21, 0, false, "TCAGG");
        auto c00_21 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_21);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_21].insert(link);
        expectedScore_[c00_21] = 1;

        link = std::make_shared<Link>();            // 00 - 31
        link->insertOccurrence(spec0, seq00, 24, false, "CTGTT");
        link->insertOccurrence(spec3, seq31, 0, false, "CTGTT");
        auto c00_31 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_31);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_31].insert(link);
        expectedScore_[c00_31] = 1;

        link = std::make_shared<Link>();            // 00 - 10 - 20
        link->insertOccurrence(spec0, seq00, 30, false, "CTTTT");
        link->insertOccurrence(spec1, seq10, 12, false, "CTTTT");
        link->insertOccurrence(spec2, seq20, 0, false, "CTTTT");
        auto c00_10_20 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_10_20);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_10_20].insert(link);
        expectedScore_[c00_10_20] = 1;

        link = std::make_shared<Link>();            // 00 - 10 - 21
        link->insertOccurrence(spec0, seq00, 36, false, "CCATC");
        link->insertOccurrence(spec1, seq10, 18, false, "CCATC");
        link->insertOccurrence(spec2, seq21, 6, false, "CCATC");
        auto c00_10_21 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_10_21);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_10_21].insert(link);
        expectedScore_[c00_10_21] = 1;

        link = std::make_shared<Link>();            // 00 - 11 - 21
        link->insertOccurrence(spec0, seq00, 42, false, "AGTTA");
        link->insertOccurrence(spec1, seq11, 6, false, "AGTTA");
        link->insertOccurrence(spec2, seq21, 12, false, "AGTTA");
        auto c00_11_21 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_11_21);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_11_21].insert(link);
        expectedScore_[c00_11_21] = 1;

        link = std::make_shared<Link>();            // 00 - 11 - 31
        link->insertOccurrence(spec0, seq00, 48, false, "TGCTG");
        link->insertOccurrence(spec1, seq11, 12, false, "TGCTG");
        link->insertOccurrence(spec3, seq31, 6, false, "TGCTG");
        auto c00_11_31 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_11_31);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_11_31].insert(link);
        expectedScore_[c00_11_31] = 1;

        link = std::make_shared<Link>();            // 00 - 32 - 40
        link->insertOccurrence(spec0, seq00, 54, false, "GGCGA");
        link->insertOccurrence(spec3, seq32, 0, false, "GGCGA");
        link->insertOccurrence(spec4, seq40, 0, false, "GGCGA");
        auto c00_32_40 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_32_40);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_32_40].insert(link);
        expectedScore_[c00_32_40] = 1;

        link = std::make_shared<Link>();            // 00 - 10 - 21 - 40
        link->insertOccurrence(spec0, seq00, 60, false, "TTCCC");
        link->insertOccurrence(spec1, seq10, 24, false, "TTCCC");
        link->insertOccurrence(spec2, seq21, 18, false, "TTCCC");
        link->insertOccurrence(spec4, seq40, 6, false, "TTCCC");
        auto c00_10_21_40 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_10_21_40);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_10_21_40].insert(link);
        expectedScore_[c00_10_21_40] = 1;

        link = std::make_shared<Link>();            // 00 - 10 - 30 - 40
        link->insertOccurrence(spec0, seq00, 66, false, "TGTAG");
        link->insertOccurrence(spec1, seq10, 30, false, "TGTAG");
        link->insertOccurrence(spec3, seq30, 0, false, "TGTAG");
        link->insertOccurrence(spec4, seq40, 12, false, "TGTAG");
        auto c00_10_30_40 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_10_30_40);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_10_30_40].insert(link);
        expectedScore_[c00_10_30_40] = 1;

        link = std::make_shared<Link>();            // 00 - 10 - 21 - 30 - 40
        link->insertOccurrence(spec0, seq00, 72, false, "CGTGT");
        link->insertOccurrence(spec1, seq10, 36, false, "CGTGT");
        link->insertOccurrence(spec2, seq21, 24, false, "CGTGT");
        link->insertOccurrence(spec3, seq30, 6, false, "CGTGT");
        link->insertOccurrence(spec4, seq40, 18, false, "CGTGT");
        auto c00_10_21_30_40 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_10_21_30_40);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_10_21_30_40].insert(link);
        expectedScore_[c00_10_21_30_40] = 1;

        link = std::make_shared<Link>();            // 00 - 11 - 21 - 31 - 40
        link->insertOccurrence(spec0, seq00, 78, false, "CGAGA");
        link->insertOccurrence(spec1, seq11, 18, false, "CGAGA");
        link->insertOccurrence(spec2, seq21, 30, false, "CGAGA");
        link->insertOccurrence(spec3, seq31, 12, false, "CGAGA");
        link->insertOccurrence(spec4, seq40, 24, false, "CGAGA");
        auto c00_11_21_31_40 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_11_21_31_40);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_11_21_31_40].insert(link);
        expectedScore_[c00_11_21_31_40] = 1;

        link = std::make_shared<Link>();            // 00 - 11 - 21 - 31 - 41
        link->insertOccurrence(spec0, seq00, 84, false, "GCCGG");
        link->insertOccurrence(spec1, seq11, 24, false, "GCCGG");
        link->insertOccurrence(spec2, seq21, 36, false, "GCCGG");
        link->insertOccurrence(spec3, seq31, 18, false, "GCCGG");
        link->insertOccurrence(spec4, seq41, 0, false, "GCCGG");
        auto c00_11_21_31_41 = std::make_shared<Cube const>(link, 100);
        expectedCubes_.insert(c00_11_21_31_41);
        expectedLinkset_.insert(link);
        expectedLinksInCubes_[c00_11_21_31_41].insert(link);
        expectedScore_[c00_11_21_31_41] = 1;

        expectedPredecessors_[c00_10_20] = {c00_10};
        expectedPredecessors_[c00_10_21] = {c00_10, c00_21};
        expectedPredecessors_[c00_11_21] = {c00_11, c00_21};
        expectedPredecessors_[c00_11_31] = {c00_11, c00_31};
        expectedPredecessors_[c00_10_30_40] = {c00_10};
        expectedPredecessors_[c00_10_21_40] = {c00_10_21};
        expectedPredecessors_[c00_10_21_30_40] = {c00_10_30_40, c00_10_21_40};

        expectedPredecessorsComplete_[c00_10_20] = {c00_10};
        expectedPredecessorsComplete_[c00_10_21] = {c00_10, c00_21};
        expectedPredecessorsComplete_[c00_11_21] = {c00_11, c00_21};
        expectedPredecessorsComplete_[c00_11_31] = {c00_11, c00_31};
        expectedPredecessorsComplete_[c00_10_30_40] = {c00_10};
        expectedPredecessorsComplete_[c00_10_21_40] = {c00_10, c00_21, c00_10_21};
        expectedPredecessorsComplete_[c00_10_21_30_40] = {c00_10, c00_21, c00_10_21, c00_10_30_40, c00_10_21_40};
    }
    auto const & expectedCubes() const { return expectedCubes_; }
    auto const & expectedLinkset() const { return expectedLinkset_; }
    auto const & expectedLinksInCubes() const { return expectedLinksInCubes_; }
    auto const & expectedPredecessors() const { return expectedPredecessors_; }
    auto const & expectedPredecessorsComplete() const { return expectedPredecessorsComplete_; }
    auto const & expectedScore() const { return expectedScore_; }

private:
    std::set<std::shared_ptr<Cube const>, CubePtrLess> expectedCubes_;
    std::unordered_set<std::shared_ptr<Link const>,
                       LinkPtrHash, LinkPtrEqual> expectedLinkset_;
    std::unordered_map<std::shared_ptr<Cube const>,
                       std::unordered_set<std::shared_ptr<Link const>,
                                          LinkPtrHash, LinkPtrEqual>,
                       CubePtrHash, CubePtrEqual> expectedLinksInCubes_;
    std::unordered_map<std::shared_ptr<Cube const>,
                       std::set<std::shared_ptr<Cube const>,
                                CubePtrLess>,
                       CubePtrHash,
                       CubePtrEqual> expectedPredecessors_;
    std::unordered_map<std::shared_ptr<Cube const>,
                       std::set<std::shared_ptr<Cube const>,
                                CubePtrLess>,
                       CubePtrHash,
                       CubePtrEqual> expectedPredecessorsComplete_;
    std::unordered_map<std::shared_ptr<Cube const>,
                       size_t,
                       CubePtrHash,
                       CubePtrEqual> expectedScore_;
};

#endif // MANUALTESTSET_H
