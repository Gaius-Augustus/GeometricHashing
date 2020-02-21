#include <array>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "catch2/catch.hpp"

#include "../DisambiguateSTLContainer.h"

TEST_CASE("Disambiguate STL Container") {
    auto array = std::array<int, 2>();
    auto map = std::map<int, int>();
    auto set = std::set<int>();
    auto unorderedMap = std::unordered_map<int, int>();
    auto unorderedSet = std::unordered_set<int>();
    auto vector = std::vector<int>();
    SECTION("Is Container") {
        REQUIRE(is_container(array));
        REQUIRE(is_container(map));
        REQUIRE(is_container(set));
        REQUIRE(is_container(unorderedMap));
        REQUIRE(is_container(unorderedSet));
        REQUIRE(is_container(vector));
    }
    SECTION("Is Associative") {
        REQUIRE(is_associative_container(map));
        REQUIRE(is_associative_container(set));
        REQUIRE(is_associative_container(unorderedMap));
        REQUIRE(is_associative_container(unorderedSet));
    }
    SECTION("Not Is Associative") {
        REQUIRE(!(is_associative_container(array)));
        REQUIRE(!(is_associative_container(vector)));
    }
    SECTION("Is Map") {
        REQUIRE(is_map_container(map));
        REQUIRE(is_map_container(unorderedMap));
    }
    SECTION("Not Is Map") {
        REQUIRE(!(is_map_container(array)));
        REQUIRE(!(is_map_container(set)));
        REQUIRE(!(is_map_container(unorderedSet)));
        REQUIRE(!(is_map_container(vector)));
    }
}
