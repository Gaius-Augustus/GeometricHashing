#include <experimental/filesystem>
#include <thread>
#include <fstream>
#include <memory>

#include "catch2/catch.hpp"
#include "../IdentifierMapping.h"
#include "../Link.h"



TEST_CASE("Link") {
    std::cout << "[INFO] -- [TEST CASE] -- Link" << std::endl;
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
        REQUIRE(link.span() == 1);
        link.extendSpanToRight(5);
        REQUIRE(link.span() == 6);
        link.extendSpanToRight(3);
        REQUIRE(link.span() == 9);
        REQUIRE(link.diagonal() == std::vector<long long>{0,0,100});

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
            auto link2 = Link(vec, 7);
        } catch (std::exception const & e) {
            thrown = true;
        }
        REQUIRE(thrown);

        vec.pop_back();
        auto link2 = Link(vec, 7);
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
        REQUIRE(link2.span() == 7);
        REQUIRE(link2.diagonal() == std::vector<long long>{0,0,100});

        REQUIRE(link2.sameDiagonal(link));
        REQUIRE(link.diagonal() == link2.diagonal());

        auto link3 = Link(std::vector<KmerOccurrence>{KmerOccurrence{0,0,10,false,"ACGT"},
                                                      KmerOccurrence{1,1,20,false,"ACGT"},
                                                      KmerOccurrence{2,2,5, false,"ACGT"},
                                                      KmerOccurrence{3,3,33,false,"ACGT"}}, 4);
        REQUIRE(link3.diagonal() == std::vector<long long>{0,10,-5,23});
    }

    SECTION("Test Link comparison") {
        Link l1;
        Link l2;
        Link l3;
        Link l4;
        Link l5;
        Link l6;
        Link l7;
        Link l8;
        Link l9;

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

        l4.insertOccurrence(0, 0, 0, false, "ACGT"); // < l1 (lower dim)
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

        l7.insertOccurrence(0, 0, 1, false, "ACGT");  // l1 < l7
        l7.insertOccurrence(1, 1, 101, false, "ACGT");
        l7.insertOccurrence(2, 1, 101, true, "ACGT");
        l7.insertOccurrence(3, 0, 501, false, "ACGT");

        l8.insertOccurrence(0, 0, 2, false, "ACGT");  // l1 < l7 < l8
        l8.insertOccurrence(1, 1, 102, false, "ACGT");
        l8.insertOccurrence(2, 1, 102, true, "ACGT");
        l8.insertOccurrence(3, 0, 502, false, "ACGT");

        l9.insertOccurrence(0, 0, 0, false, "AGCT");  // == l1, == l5
        l9.insertOccurrence(1, 1, 100, false, "AGCT");
        l9.insertOccurrence(2, 1, 100, true, "AGCT");
        l9.insertOccurrence(3, 0, 500, false, "AGCT");

        REQUIRE(l1 == l2);
        REQUIRE(!(l1 < l2));
        REQUIRE(!(l2 < l1));
        REQUIRE(!(l1 < l3));
        REQUIRE(l3 < l1);
        REQUIRE(!(l1 == l3));
        REQUIRE(l4 < l1);
        REQUIRE(!(l1 < l4));
        REQUIRE(!(l1 == l4));
        REQUIRE(l1 == l5);
        REQUIRE(!(l5 < l1));
        REQUIRE(!(l1 < l5));
        REQUIRE(!(l1 < l6));
        REQUIRE(l6 < l1);
        REQUIRE(!(l1 == l6));
        REQUIRE(l5 < l7);
        REQUIRE(!(l7 < l5));
        REQUIRE(!(l5 == l7));
        REQUIRE(l7 < l8);
        REQUIRE(!(l8 < l7));
        REQUIRE(!(l7 == l8));
        REQUIRE(l1 == l9);
        REQUIRE(!(l9 < l1));
        REQUIRE(!(l1 < l9));
    }
}



TEST_CASE("Link Pair") {
    std::cout << "[INFO] -- [TEST CASE] -- Link Pair" << std::endl;
    KmerOccurrence occ1_1(0,1,100,false,"AAAA");
    KmerOccurrence occ1_2(1,2,50,false,"AAAA");
    Link dist1{std::vector<KmerOccurrence>{occ1_2, occ1_1}, 4};   // should store occ1 as first
    REQUIRE(dist1.first() == occ1_1);
    REQUIRE(dist1.second() == occ1_2);
    REQUIRE(dist1.diagonal().at(1) == -50);

    KmerOccurrence occ2_1(0,0,100,false,"AAAA");  // smaller seq
    KmerOccurrence occ2_2(1,2,50,false,"AAAA");
    Link dist2{std::vector<KmerOccurrence>{occ2_1, occ2_2}, 4};
    REQUIRE(dist2.first() == occ2_1);
    REQUIRE(dist2.second() == occ2_2);
    REQUIRE(dist2.diagonal().at(1) == -50);
    REQUIRE(dist2.first() < dist1.first());
    REQUIRE(dist2.second() == dist1.second());
    REQUIRE(dist2 < dist1);
    REQUIRE(!(dist1 < dist2));
    REQUIRE(!(dist1.sameSequences(dist2)));
    REQUIRE(!(dist1.sameDiagonal(dist2)));

    KmerOccurrence occ3_1(0,1,90,false,"AAAA");  // bigger distance, smaller pos
    KmerOccurrence occ3_2(1,2,50,false,"AAAA");
    Link dist3{std::vector<KmerOccurrence>{occ3_1, occ3_2}, 4};
    REQUIRE(dist3.diagonal().at(1) == -40);
    REQUIRE(dist3.first() < dist1.first());
    REQUIRE(dist3.second() == dist1.second());
    REQUIRE(dist1 < dist3); // because distance gets evaluated before positions!
    REQUIRE(dist1.sameSequences(dist3));
    REQUIRE(!(dist1.sameDiagonal(dist3)));

    KmerOccurrence occ4_1(0,1,100,false,"AAAA");
    KmerOccurrence occ4_2(1,2,40,false,"AAAA");   // smaller distance, smaller pos
    Link dist4{std::vector<KmerOccurrence>{occ4_1, occ4_2}, 4};
    REQUIRE(dist4.diagonal().at(1) == -60);
    REQUIRE(dist4.first() == dist1.first());
    REQUIRE(dist4.second() < dist1.second());
    REQUIRE(dist4 < dist1);
    REQUIRE(dist1.sameSequences(dist4));
    REQUIRE(!(dist1.sameDiagonal(dist4)));

    KmerOccurrence occ5_1(0,1,90,false,"AAAA");  // same distance, smaller pos
    KmerOccurrence occ5_2(1,2,40,false,"AAAA");
    Link dist5{std::vector<KmerOccurrence>{occ5_1, occ5_2}, 4};
    REQUIRE(dist5.first() < dist1.first());
    REQUIRE(dist5.second() < dist1.second());
    REQUIRE(dist1.sameSequences(dist5));
    REQUIRE(dist1.sameDiagonal(dist5));
    REQUIRE(dist5.sameDiagonal(dist1));

    SECTION("Test Diagonal Sorting") {
        Link diag1_1{std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{1,1,10,false,"AAAA"}}, 4};
        Link diag1_2{std::vector<KmerOccurrence>{KmerOccurrence{0,0,1,false,"AAAA"}, KmerOccurrence{1,1,11,false,"AAAA"}}, 4};
        Link diag1_3{std::vector<KmerOccurrence>{KmerOccurrence{0,0,2,false,"AAAA"}, KmerOccurrence{1,1,12,false,"AAAA"}}, 4};
        Link diag2_1{std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{1,1,0,false,"AAAA"}}, 4};
        Link diag2_2{std::vector<KmerOccurrence>{KmerOccurrence{0,0,1,false,"AAAA"}, KmerOccurrence{1,1,1,false,"AAAA"}}, 4};
        Link diag2_3{std::vector<KmerOccurrence>{KmerOccurrence{0,0,2,false,"AAAA"}, KmerOccurrence{1,1,2,false,"AAAA"}}, 4};
        Link diag3_1{std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{2,2,10,false,"AAAA"}}, 4};
        Link diag3_2{std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{2,2,11,false,"AAAA"}}, 4};
        Link diag3_3{std::vector<KmerOccurrence>{KmerOccurrence{0,0,0,false,"AAAA"}, KmerOccurrence{2,2,12,false,"AAAA"}}, 4};
        std::vector<Link> unsorted{diag3_3, diag2_1, diag1_2,
                                   diag1_1, diag2_3, diag1_3,
                                   diag2_2, diag3_2, diag3_1};
        std::vector<Link> sorted{diag2_1, diag2_2, diag2_3,
                                 diag1_1, diag1_2, diag1_3,
                                 diag3_1, diag3_2, diag3_3};
        std::sort(unsorted.begin(), unsorted.end());
        REQUIRE(unsorted == sorted);
    }
}
