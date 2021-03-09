#include <bitset>
#include <vector>

#include "catch2/catch.hpp"
#include "../SpacedSeedMask.h"
#include "../SpacedSeedMaskCollection.h"

TEST_CASE("SpacedSeedMask") {
    std::cout << "[INFO] -- [TEST CASE] -- SpacedSeedMask" << std::endl;
    SpacedSeedMask mask(42, 120);
    REQUIRE(mask.span() == 120);
    REQUIRE(mask.weight() == 42);

    SECTION("Test rng") {
        SpacedSeedMask mask2(42, 120);
        REQUIRE(mask.mask() != mask2.mask());
    }

    SECTION("Test throw") {
        auto thrown = false;
        try {
            SpacedSeedMask mask3(121, 120);
        } catch (std::runtime_error const & e) {
            thrown = true;
        }
        REQUIRE(thrown);
    }

    SECTION("Test set bits vector") {
        SpacedSeedMask mask4(42, 64);
        auto setBitsVector = mask4.getSetPositions();
        std::bitset<64> stlBitset;
        for (auto&& pos : setBitsVector) {
            stlBitset.set(pos);
        }
        REQUIRE(stlBitset.to_ulong() == mask4.mask().to_ulong());
    }

    SECTION("Test String Construction") {
        SpacedSeedMask mask("0101001100");
        REQUIRE(mask.span() == 10);
        REQUIRE(mask.weight() == 4);
        REQUIRE(mask.getSetPositions() == std::vector<size_t>{1,3,6,7});
    }

    SECTION("Test String Method") {
        std::string maskStr{"0101001100"};
        SpacedSeedMask mask(maskStr);
        REQUIRE(mask.string() == maskStr);
    }
}



TEST_CASE("SpacedSeedMaskCollection") {
    std::cout << "[INFO] -- [TEST CASE] -- SpacedSeedMaskCollection" << std::endl;
    auto masks = SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3),
                                          SpacedSeedMaskCollection::Span(6),
                                          SpacedSeedMaskCollection::SeedSetSize(20));    // all possible masks
    REQUIRE(masks.maxSpan() == 6);
    REQUIRE(masks.weight() == 3);
    REQUIRE(masks.size() == 20);
    for (auto it = masks.masks().begin(); it != masks.masks().end(); ++it) {
        auto it2 = it;
        ++it2;
        for (; it2 != masks.masks().end(); ++it2) {
            REQUIRE(!(*it == *it2));
        }
    }

    SECTION("Test String Construction") {
        auto masks = SpacedSeedMaskCollection({"0101001100", "0010010101"});
        REQUIRE(masks.maxSpan() == 10);
        REQUIRE(masks.weight() == 4);
        REQUIRE(masks.size() == 2);
        auto masks2 = SpacedSeedMaskCollection({"0101001100", "00100101010"});
        REQUIRE(masks2.maxSpan() == 11);
        REQUIRE(masks2.span(0) == 10);
        REQUIRE(masks2.span(1) == 11);
        REQUIRE(masks2.weight() == 4);
        REQUIRE(masks2.size() == 2);
    }
}

TEST_CASE("Optimal Seeds Single") {
    std::cout << "[INFO] -- [TEST CASE] -- Optimal Seeds Single" << std::endl;
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3)).masks().at(0) == SpacedSeedMask("1101"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3)).maxSpan() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3)).weight() == 3);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4)).masks().at(0) == SpacedSeedMask("1100101"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4)).maxSpan() == 7);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4)).weight() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5)).masks().at(0) == SpacedSeedMask("10010111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5)).maxSpan() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5)).weight() == 5);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6)).masks().at(0) == SpacedSeedMask("1100101011"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6)).maxSpan() == 10);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6)).weight() == 6);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7)).masks().at(0) == SpacedSeedMask("11101010011"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7)).maxSpan() == 11);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7)).weight() == 7);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8)).masks().at(0) == SpacedSeedMask("11001001010111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8)).maxSpan() == 14);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8)).weight() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9)).masks().at(0) == SpacedSeedMask("111000101011011"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9)).maxSpan() == 15);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9)).weight() == 9);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10)).masks().at(0) == SpacedSeedMask("1101100011010111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10)).maxSpan() == 16);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10)).weight() == 10);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11)).masks().at(0) == SpacedSeedMask("111001011001010111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11)).maxSpan() == 18);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11)).weight() == 11);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12)).masks().at(0) == SpacedSeedMask("111011001011010111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12)).maxSpan() == 18);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12)).weight() == 12);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13)).masks().at(0) == SpacedSeedMask("11101110010110010111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13)).maxSpan() == 20);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13)).weight() == 13);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14)).masks().at(0) == SpacedSeedMask("111101001101001110111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14)).maxSpan() == 21);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14)).weight() == 14);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15)).masks().at(0) == SpacedSeedMask("11110010101011001101111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15)).maxSpan() == 23);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15)).weight() == 15);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16)).masks().at(0) == SpacedSeedMask("111100110101011001101111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16)).maxSpan() == 24);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16)).weight() == 16);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17)).masks().at(0) == SpacedSeedMask("1111010110011001101011111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17)).maxSpan() == 25);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17)).weight() == 17);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18)).masks().at(0) == SpacedSeedMask("111101100101101010110011111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18)).maxSpan() == 27);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18)).weight() == 18);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19)).masks().at(0) == SpacedSeedMask("11111011001110101101011111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19)).maxSpan() == 26);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19)).weight() == 19);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20)).masks().at(0) == SpacedSeedMask("11110110100111010011101110111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20)).maxSpan() == 29);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20)).weight() == 20);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21)).masks().at(0) == SpacedSeedMask("111110101100111001011101101111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21)).maxSpan() == 30);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21)).weight() == 21);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22)).masks().at(0) == SpacedSeedMask("111110111001110101011011011111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22)).maxSpan() == 30);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22)).weight() == 22);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23)).masks().at(0) == SpacedSeedMask("111111001101011010111001101011111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23)).maxSpan() == 33);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23)).weight() == 23);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24)).masks().at(0) == SpacedSeedMask("1111110101011100111001011011011111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24)).maxSpan() == 34);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24)).weight() == 24);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(25)).masks().at(0) == SpacedSeedMask("111110111011100111010110110111111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(25)).maxSpan() == 33);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(25)).weight() == 25);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(25)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(26)).masks().at(0) == SpacedSeedMask("11111101011100111011011010111011111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(26)).maxSpan() == 35);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(26)).weight() == 26);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(26)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(27)).masks().at(0) == SpacedSeedMask("111111011010111010111001110110111111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(27)).maxSpan() == 36);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(27)).weight() == 27);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(27)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(28)).masks().at(0) == SpacedSeedMask("1111110110101110110011110101110111111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(28)).maxSpan() == 37);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(28)).weight() == 28);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(28)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(29)).masks().at(0) == SpacedSeedMask("11111011101101110011110110101110111111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(29)).maxSpan() == 38);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(29)).weight() == 29);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(29)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(30)).masks().at(0) == SpacedSeedMask("11111011110110111010111011110110111111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(30)).maxSpan() == 38);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(30)).weight() == 30);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(30)).size() == 1);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(31)).masks().at(0) == SpacedSeedMask("111111011110111010111011011110110111111"));
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(31)).maxSpan() == 39);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(31)).weight() == 31);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(31)).size() == 1);
}

TEST_CASE("Optimal Seeds Set") {
    std::cout << "[INFO] -- [TEST CASE] -- Optimal Seeds Set" << std::endl;
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("101001"), SpacedSeedMask("1000000000011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 13);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 3);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1100101"), SpacedSeedMask("11000000101") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 11);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("10110011"), SpacedSeedMask("101001000011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 12);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 5);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1100101011"), SpacedSeedMask("110001000001101") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 15);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 6);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11101001101"), SpacedSeedMask("111000100001011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 15);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 7);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("110110100111"), SpacedSeedMask("10110001000100111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 17);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110010111"), SpacedSeedMask("110010100000010100111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 21);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 9);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111010010110111"), SpacedSeedMask("11101001000100010111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 20);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 10);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110100110111"), SpacedSeedMask("111010001100010010111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 21);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 11);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111010110100110111"), SpacedSeedMask("11110010001000101001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 23);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 12);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1111011001011010111"), SpacedSeedMask("1110100100010011000101111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 25);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 13);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111010110110011001111"), SpacedSeedMask("1110110000101000100101001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 28);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 14);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111101100110101011111"), SpacedSeedMask("1110110010100011000010110111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 28);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 15);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1111011011100101110111"), SpacedSeedMask("111101010001011000010011001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 30);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 16);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111101011011001101011111"), SpacedSeedMask("111110101000100001010010011001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 33);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 17);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1111101011010111001101111"), SpacedSeedMask("1111101001100010001010000110110111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 34);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 18);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11111010110111001011101111"), SpacedSeedMask("11110111000100100010100010101101111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 35);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 19);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111011110110110011101011111"), SpacedSeedMask("11111011000101001000110001011011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 35);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 20);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11110110101011011001111011111"), SpacedSeedMask("11111011001001100011000010010101011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 38);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 21);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11111101101101011100111011111"), SpacedSeedMask("111110101100010010100010100100111011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 39);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 22);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111110101111011100110110111111"), SpacedSeedMask("11111011010001010100100111000110111111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 38);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 23);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(2)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11111011101101011100111011011111"), SpacedSeedMask("111110111000101100101001011000110111111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(2)).maxSpan() == 39);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(2)).weight() == 24);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(2)).size() == 2);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("10011"), SpacedSeedMask("101000001"), SpacedSeedMask("10000000011"), SpacedSeedMask("1000010000001") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 13);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 3);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("101011"), SpacedSeedMask("11000101"), SpacedSeedMask("1010000011"), SpacedSeedMask("10010000000011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 14);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11101001"), SpacedSeedMask("1001100011"), SpacedSeedMask("1010001011"), SpacedSeedMask("11000000100101") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 14);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 5);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111010011"), SpacedSeedMask("110001001011"), SpacedSeedMask("11001000001011"), SpacedSeedMask("101010000001000011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 18);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 6);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1110101011"), SpacedSeedMask("110100100111"), SpacedSeedMask("1100010001001011"), SpacedSeedMask("110010000001010011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 18);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 7);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("110101011011"), SpacedSeedMask("111001000101011"), SpacedSeedMask("111000101000010011"), SpacedSeedMask("110100100000010000111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 21);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110100111"), SpacedSeedMask("11010000110010111"), SpacedSeedMask("11100010010000101011"), SpacedSeedMask("11010001000001000100111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 23);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 9);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11101010110111"), SpacedSeedMask("110110010100010111"), SpacedSeedMask("11101001000100001010011"), SpacedSeedMask("11010010000100001001000111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 26);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 10);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111011011010111"), SpacedSeedMask("111010100011001111"), SpacedSeedMask("1110100010010100100111"), SpacedSeedMask("11110001001000000100010111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 26);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 11);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110110101111"), SpacedSeedMask("1111010100011010111"), SpacedSeedMask("111100010010000100110111"), SpacedSeedMask("11101000100010000010010001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 29);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 12);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110101001101111"), SpacedSeedMask("111010100010110000111011"), SpacedSeedMask("111010011100001001010111"), SpacedSeedMask("111100100010000001000100101111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 30);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 13);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111101101011101111"), SpacedSeedMask("11101100110010010101111"), SpacedSeedMask("111100011001010001001001111"), SpacedSeedMask("1110100101100000001000110010111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 31);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 14);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11101110110110101111"), SpacedSeedMask("1111101001010001100110111"), SpacedSeedMask("111101000100010010010100011111"), SpacedSeedMask("111100011001000100001000101001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 33);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 15);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111110110011110101111"), SpacedSeedMask("1111001001010001100101011111"), SpacedSeedMask("111101010100001101001001100111"), SpacedSeedMask("11101100010010001000010100001101111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 35);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 16);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1111010111011101101111"), SpacedSeedMask("11110110001101001010001011111"), SpacedSeedMask("111101000101001001000101100011111"), SpacedSeedMask("11110010110000010001000100100001110111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 38);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 17);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11101111010110110111111"), SpacedSeedMask("11110110011000110101010011111"), SpacedSeedMask("111110010010101001001001100011111"), SpacedSeedMask("111110001100100100000101000100001011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 39);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 18);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11110110111010111001101111"), SpacedSeedMask("1111001010011100101010011011111"), SpacedSeedMask("11111001100100001001101001010101111"), SpacedSeedMask("1111010100100010100010000011001001011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 40);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 19);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11111011010111011101101111"), SpacedSeedMask("111110110101001100001100101011111"), SpacedSeedMask("111101100010100001100010010010101011111"), SpacedSeedMask("1111100101001010001000000101100001011001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 43);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 20);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111110110111101011101101111"), SpacedSeedMask("111110110001100101011010100111111"), SpacedSeedMask("11111010100011010010000110010011011111"), SpacedSeedMask("1111100110100010000010100001100101001011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 43);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 21);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111110111101101110111011111"), SpacedSeedMask("111111011000111000100101001101011111"), SpacedSeedMask("11101011101100001011010001011001101111"), SpacedSeedMask("11111100101000100100100000110010100011011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 44);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 22);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11111011110111011101101101111"), SpacedSeedMask("11111011011000101001110011010111111"), SpacedSeedMask("11101101100001111101001000110101011111"), SpacedSeedMask("1111011100011010100000011001001001011011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 43);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 23);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(4)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111110111101011111011011110111"), SpacedSeedMask("1111101101011010110001001110010111111"), SpacedSeedMask("111101110011001000101110100101110101111"), SpacedSeedMask("111111001100100101000100000110010101001111111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(4)).maxSpan() == 45);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(4)).weight() == 24);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(4)).size() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1011"), SpacedSeedMask("1010000001"), SpacedSeedMask("1000000000011"), SpacedSeedMask("10010000000001"), SpacedSeedMask("100010000000001"), SpacedSeedMask("100001000000001"), SpacedSeedMask("100000100000001"), SpacedSeedMask("100000001000001") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 15);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 3);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(3), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("10111"), SpacedSeedMask("10001000011"), SpacedSeedMask("10010000000101"), SpacedSeedMask("100000000100101"), SpacedSeedMask("1000001000000101"), SpacedSeedMask("1001000010000001"), SpacedSeedMask("1000100000000011"), SpacedSeedMask("1000000010000011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 16);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 4);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(4), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("110111"), SpacedSeedMask("110001011"), SpacedSeedMask("101001000011"), SpacedSeedMask("1010001001001"), SpacedSeedMask("10100100000011"), SpacedSeedMask("100100010000101"), SpacedSeedMask("110000100010001"), SpacedSeedMask("1010000010000000011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 19);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 5);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(5), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11010111"), SpacedSeedMask("11100100101"), SpacedSeedMask("110010001011"), SpacedSeedMask("1011000010011"), SpacedSeedMask("11000010100011"), SpacedSeedMask("1010000100010011"), SpacedSeedMask("1100100000000010101"), SpacedSeedMask("11000010000000001101") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 20);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 6);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(6), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("110101111"), SpacedSeedMask("110101000111"), SpacedSeedMask("11100100010011"), SpacedSeedMask("1100101000010101"), SpacedSeedMask("11010001000001011"), SpacedSeedMask("11000100100000100011"), SpacedSeedMask("110100000000100001011"), SpacedSeedMask("1010010000010000010011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 22);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 7);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(7), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("110100111011"), SpacedSeedMask("11101010001011"), SpacedSeedMask("110101001000111"), SpacedSeedMask("110010010010000111"), SpacedSeedMask("1101000001000101011"), SpacedSeedMask("1100101000010010011"), SpacedSeedMask("110101000010000010011"), SpacedSeedMask("11100010000000100001011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 23);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(8), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11110110111"), SpacedSeedMask("11101010001111"), SpacedSeedMask("1101100010010111"), SpacedSeedMask("11100001010001010011"), SpacedSeedMask("110010010000100011011"), SpacedSeedMask("101011000001001000111"), SpacedSeedMask("1100101000100001001011"), SpacedSeedMask("111000010010000001000111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 24);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 9);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(9), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111110110111"), SpacedSeedMask("11110001010011011"), SpacedSeedMask("110101100001010111"), SpacedSeedMask("11010100010011000111"), SpacedSeedMask("1101000010100101000111"), SpacedSeedMask("11100100100000100101011"), SpacedSeedMask("11100010001000000100100111"), SpacedSeedMask("11101000000100010001001011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 26);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 10);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(10), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11101011011111"), SpacedSeedMask("11101001101010111"), SpacedSeedMask("111001010110000010111"), SpacedSeedMask("111010100000110100111"), SpacedSeedMask("1101100010010001101011"), SpacedSeedMask("111000101001001000100111"), SpacedSeedMask("110110001100010001001011"), SpacedSeedMask("111100010000010000100010111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 27);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 11);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(11), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111011101101111"), SpacedSeedMask("11101100001100110111"), SpacedSeedMask("111100011010001010111"), SpacedSeedMask("111010010010100100001111"), SpacedSeedMask("111001010101000011001011"), SpacedSeedMask("11101000110000010010100111"), SpacedSeedMask("1101101000001000100001010111"), SpacedSeedMask("11100010010000100010010010111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 29);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 12);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(12), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1111011011101111"), SpacedSeedMask("11101011100010110111"), SpacedSeedMask("111101000101101001111"), SpacedSeedMask("1101100011010000101010111"), SpacedSeedMask("110101110000100010000110111"), SpacedSeedMask("111010010010001010001001111"), SpacedSeedMask("1110110001000011000100100111"), SpacedSeedMask("111001010010010000100100010111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 30);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 13);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(13), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111101101011101111"), SpacedSeedMask("111101011100010110111"), SpacedSeedMask("111010100001011001101111"), SpacedSeedMask("1110110001101000101010111"), SpacedSeedMask("11101001100011001100011011"), SpacedSeedMask("111001011001000100101001111"), SpacedSeedMask("11101010010001000001000100101111"), SpacedSeedMask("11101101000010010000010000111011") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 32);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 14);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(14), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111111011101101111"), SpacedSeedMask("1111101001011100110111"), SpacedSeedMask("1101010111000100011101111"), SpacedSeedMask("11110011000110100011010111"), SpacedSeedMask("1111001001000101101000101111"), SpacedSeedMask("11101101001010001010010010111"), SpacedSeedMask("11110100010100001000100110001111"), SpacedSeedMask("111101010000100100000010010110111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 33);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 15);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(15), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1111110111101110111"), SpacedSeedMask("111011100011100101101111"), SpacedSeedMask("1111010011000101110110111"), SpacedSeedMask("1111010010011010010001011111"), SpacedSeedMask("11110011001010000101100001011011"), SpacedSeedMask("110111001000010100100010110001111"), SpacedSeedMask("11110010101000100000001100101001111"), SpacedSeedMask("11110110010010000010100001001100111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 35);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 16);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(16), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11111101111011101111"), SpacedSeedMask("1111101011001011100101111"), SpacedSeedMask("111010110001101001101011111"), SpacedSeedMask("111101001101000010101001110111"), SpacedSeedMask("11011101000010110010000110101111"), SpacedSeedMask("1111100100101000100011010001001111"), SpacedSeedMask("111100110010000100101000110010001111"), SpacedSeedMask("111100101000110000100010000110110111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 36);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 17);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(17), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1111101110110111101111"), SpacedSeedMask("111110110101100101100110111"), SpacedSeedMask("1110110100111010000110001011111"), SpacedSeedMask("111101010001100001010011110010111"), SpacedSeedMask("111101010100000110110001001101111"), SpacedSeedMask("1110011011001010000001010101001011011"), SpacedSeedMask("1111010000110010100100100010010101111"), SpacedSeedMask("1111000110001001000110010100110001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 37);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 18);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(18), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111110110110111011110111"), SpacedSeedMask("11111001110010101100101011111"), SpacedSeedMask("111011101010011000011110110111"), SpacedSeedMask("111010111001000101110010010011111"), SpacedSeedMask("11110100100110100010001010110011111"), SpacedSeedMask("111101011000100101000101000010001110111"), SpacedSeedMask("111100100111000000101001001101000110111"), SpacedSeedMask("111110100000101001010000101001011001111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 39);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 19);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(19), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111101111011011111011111"), SpacedSeedMask("11110110010111011001010111111"), SpacedSeedMask("111110011101010010111001101111"), SpacedSeedMask("11111100100010110001110100100101111"), SpacedSeedMask("11110101011100100100000110011011111"), SpacedSeedMask("1111010110001100000011001001010010011111"), SpacedSeedMask("1111100100101001100010000110001010101111"), SpacedSeedMask("1110110011010000010100110100001010110111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 40);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 20);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(20), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11111011101110110111101111"), SpacedSeedMask("1111101001101001100111010111111"), SpacedSeedMask("1111011001010101101011100100011111"), SpacedSeedMask("111010101100111110000010100110110111"), SpacedSeedMask("111110111001010000010100100011011010111"), SpacedSeedMask("111101101001100011001000000011100010101111"), SpacedSeedMask("111110001011010000001100010100001101110111"), SpacedSeedMask("111101001100000101100001101010001001011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 42);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 21);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(21), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("111011011101111101111011111"), SpacedSeedMask("11111001101111001001101010111111"), SpacedSeedMask("1111101101000100111010110110011111"), SpacedSeedMask("1111001101010111000110001001011011111"), SpacedSeedMask("1110111101100000010101100101010010110111"), SpacedSeedMask("1111100101010110001001000010100011000111111"), SpacedSeedMask("1111001011000010110100010000111000101110111"), SpacedSeedMask("1110111010100001100000011011000110010110111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 43);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 22);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(22), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("1111011110111101110110111111"), SpacedSeedMask("11111011001101010111110011101111"), SpacedSeedMask("11110100111100111110010100101011111"), SpacedSeedMask("111011110000111011000101011101011111"), SpacedSeedMask("11111101010100100001101110010011011111"), SpacedSeedMask("111110100110011010100000100010110100111111"), SpacedSeedMask("111110101100101001000010011000001100111001111"), SpacedSeedMask("111111010011000010001011000011000100101011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 45);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 23);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(23), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(8)).masks() == std::vector<SpacedSeedMask>{ SpacedSeedMask("11111101101110111011110111111"), SpacedSeedMask("1111100111010110100101110110111111"), SpacedSeedMask("11110110101101101111100100011101111"), SpacedSeedMask("11111001111001000110010110101101011111"), SpacedSeedMask("11110111001100001111010001010010010111111"), SpacedSeedMask("11110100110001110010101011100001110101111"), SpacedSeedMask("11111101010010011000010010001011100011101111"), SpacedSeedMask("1111101001101010000001100011000100110110011111") });
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(8)).maxSpan() == 46);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(8)).weight() == 24);
    REQUIRE(SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(8)).size() == 8);

    auto thrown = false;
    try {
        SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(2));
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(32));
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
    thrown = false;
    try {
        SpacedSeedMaskCollection(SpacedSeedMaskCollection::Weight(24), SpacedSeedMaskCollection::SeedSetSize(3));
    } catch (std::exception const & e) {
        thrown = true;
    }
    REQUIRE(thrown);
}
