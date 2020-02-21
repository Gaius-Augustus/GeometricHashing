#include <memory>
#include <unordered_set>

#include "catch2/catch.hpp"
#include "../TwoBitKmer.h"

TEST_CASE("TwoBitKmer") {
    // less than 29 bases, 29 bases, more than 29 but less than 61 bases, more than 61 bases
    for (auto&& kmerLiteral : {"ACGT", "AGTGTGTTATATGAGCATATAATGCATGT",
                               "AGTGTGTTATATGAGCATATAATGCATGTGTGTATTCTCGAGGGCTGAGG",
                               "AGTGTGTTATATGAGCATATAATGCATGTGTGTATTCTCGAGGGCTGAGAGTGTGTTATATGAGC"}) {
        std::string kmerString{kmerLiteral};

        SECTION("Test basic class functionality") {
            if (kmerString.size() <= 29) {
                auto kmer = std::make_shared<TwoBitKmer<TwoBitKmerDataShort>>(kmerString);
                auto kmer2 = std::make_shared<TwoBitKmer<TwoBitKmerDataShort>>(kmerLiteral);

                REQUIRE(kmer->length() == kmerString.length());
                REQUIRE(kmer->toString() == kmerString);
                REQUIRE(*kmer == kmerString);
                REQUIRE(*kmer == *kmer2);
                for (size_t i = 0; i < kmerString.length(); ++i) {
                    REQUIRE(kmer->base(i) == kmerString.at(i));
                }
            } else if (kmerString.size() <= 61) {
                auto kmer = std::make_shared<TwoBitKmer<TwoBitKmerDataMedium>>(kmerString);
                auto kmer2 = std::make_shared<TwoBitKmer<TwoBitKmerDataMedium>>(kmerLiteral);

                REQUIRE(kmer->length() == kmerString.length());
                REQUIRE(kmer->toString() == kmerString);
                REQUIRE(*kmer == kmerString);
                REQUIRE(*kmer == *kmer2);
                for (size_t i = 0; i < kmerString.length(); ++i) {
                    REQUIRE(kmer->base(i) == kmerString.at(i));
                }
            } else {
                auto kmer = std::make_shared<TwoBitKmer<TwoBitKmerDataLong>>(kmerString);
                auto kmer2 = std::make_shared<TwoBitKmer<TwoBitKmerDataLong>>(kmerLiteral);

                REQUIRE(kmer->length() == kmerString.length());
                REQUIRE(kmer->toString() == kmerString);
                REQUIRE(*kmer == kmerString);
                REQUIRE(*kmer == *kmer2);
                for (size_t i = 0; i < kmerString.length(); ++i) {
                    REQUIRE(kmer->base(i) == kmerString.at(i));
                }
            }
        }

        SECTION("Test reverse complement stuff") {
            if (kmerString.size() <= 29) {
                TwoBitKmer<TwoBitKmerDataShort> kmer(kmerString);
                auto kmerRC = kmer.reverseComplement();
                TwoBitKmer<TwoBitKmerDataShort> kmerRC2{reverseComplement(kmerString)};
                REQUIRE(kmerRC == kmerRC2);

                TwoBitKmerRCIncludingHash<TwoBitKmerDataShort> hasher;
                REQUIRE(hasher(kmer) == hasher(kmerRC));
            } else if (kmerString.size() <= 61) {
                TwoBitKmer<TwoBitKmerDataMedium> kmer(kmerString);
                auto kmerRC = kmer.reverseComplement();
                TwoBitKmer<TwoBitKmerDataMedium> kmerRC2{reverseComplement(kmerString)};
                REQUIRE(kmerRC == kmerRC2);

                TwoBitKmerRCIncludingHash<TwoBitKmerDataMedium> hasher;
                REQUIRE(hasher(kmer) == hasher(kmerRC));
            } else {
                TwoBitKmer<TwoBitKmerDataLong> kmer(kmerString);
                auto kmerRC = kmer.reverseComplement();
                TwoBitKmer<TwoBitKmerDataLong> kmerRC2{reverseComplement(kmerString)};
                REQUIRE(kmerRC == kmerRC2);

                TwoBitKmerRCIncludingHash<TwoBitKmerDataLong> hasher;
                REQUIRE(hasher(kmer) == hasher(kmerRC));
            }
        }
    }

    SECTION("Test comparisons") {
        TwoBitKmer<TwoBitKmerDataShort> kmer{"AGGTCTAGACCT"};    // palindrome
        TwoBitKmer<TwoBitKmerDataShort> kmer2{"TGGTCTAGACCT"};   // reverse complement: AGGTCTAGACCA which is less
        TwoBitKmerRCIncludingLess<TwoBitKmerDataShort> compareLess;
        REQUIRE(kmer == kmer.reverseComplement());
        REQUIRE(!(kmer2 == kmer2.reverseComplement()));
        REQUIRE(compareLess(kmer2, kmer));
    }

    SECTION("Test RCIncluding Set (short)") {
        TwoBitKmer<TwoBitKmerDataShort> kmer{"AGTGTGTTATATGAGCATATAATGCATGT"};
        auto kmerRC = kmer.reverseComplement();

        std::unordered_set<TwoBitKmer<TwoBitKmerDataShort>,
                           TwoBitKmerRCIncludingHash<TwoBitKmerDataShort>,
                           TwoBitKmerRCIncludingEqual<TwoBitKmerDataShort>> kmerSet;
        kmerSet.emplace(kmer);
        REQUIRE(kmerSet.find(kmerRC) != kmerSet.end());
        auto insertPair = kmerSet.emplace(kmerRC);
        REQUIRE(*(kmerSet.begin()) == kmer);
        REQUIRE(insertPair.second == false);    // no insertion took place
    }

    SECTION("Test RCIncluding Set (medium)") {
        TwoBitKmer<TwoBitKmerDataMedium> kmer{"AGTGTGTTATATGAGCATATAATGCATGTGTGTATTCTCGAGGGCTGAGG"};
        auto kmerRC = kmer.reverseComplement();

        std::unordered_set<TwoBitKmer<TwoBitKmerDataMedium>,
                           TwoBitKmerRCIncludingHash<TwoBitKmerDataMedium>,
                           TwoBitKmerRCIncludingEqual<TwoBitKmerDataMedium>> kmerSet;
        kmerSet.emplace(kmer);
        REQUIRE(kmerSet.find(kmerRC) != kmerSet.end());
        auto insertPair = kmerSet.emplace(kmerRC);
        REQUIRE(*(kmerSet.begin()) == kmer);
        REQUIRE(insertPair.second == false);    // no insertion took place
    }

    SECTION("Test RCIncluding Set (long)") {
        TwoBitKmer<TwoBitKmerDataLong> kmer{"AGTGTGTTATATGAGCATATAATGCATGTGTGTATTCTCGAGGGCTGAGGAGTGTGTTATATGAGCATATAATGCATGTGTGTATTCTCGAGGGCTGAGG"};
        auto kmerRC = kmer.reverseComplement();

        std::unordered_set<TwoBitKmer<TwoBitKmerDataLong>,
                           TwoBitKmerRCIncludingHash<TwoBitKmerDataLong>,
                           TwoBitKmerRCIncludingEqual<TwoBitKmerDataLong>> kmerSet;
        kmerSet.emplace(kmer);
        REQUIRE(kmerSet.find(kmerRC) != kmerSet.end());
        auto insertPair = kmerSet.emplace(kmerRC);
        REQUIRE(*(kmerSet.begin()) == kmer);
        REQUIRE(insertPair.second == false);    // no insertion took place
    }
}
