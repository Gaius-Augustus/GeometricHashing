#include <array>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "catch2/catch.hpp"
#include "../ContainerChunks.h"

using namespace mabl3;



TEST_CASE("Test getContainerChunk()") {
    std::cout << "[INFO] -- [TEST CASE] -- Test getContainerChunk()" << std::endl;
    SECTION("Vector") {
        std::vector<int> container({0,1,2,3,4,5,6,7,8,9});
        SECTION("More Threads Than Elements") {
            REQUIRE(getContainerChunk(container, 0, 11) == std::vector<int>{0});
            REQUIRE(getContainerChunk(container, 1, 11) == std::vector<int>{1});
            REQUIRE(getContainerChunk(container, 2, 11) == std::vector<int>{2});
            REQUIRE(getContainerChunk(container, 3, 11) == std::vector<int>{3});
            REQUIRE(getContainerChunk(container, 4, 11) == std::vector<int>{4});
            REQUIRE(getContainerChunk(container, 5, 11) == std::vector<int>{5});
            REQUIRE(getContainerChunk(container, 6, 11) == std::vector<int>{6});
            REQUIRE(getContainerChunk(container, 7, 11) == std::vector<int>{7});
            REQUIRE(getContainerChunk(container, 8, 11) == std::vector<int>{8});
            REQUIRE(getContainerChunk(container, 9, 11) == std::vector<int>{9});
            REQUIRE(getContainerChunk(container, 10, 11) == std::vector<int>{});
            auto thrown = false;
            try {
                getContainerChunk(container, 11, 11);
            } catch (std::exception const & e) {
                thrown = true;
            }
            REQUIRE(thrown);
        }
        SECTION("Even Chunks") {
            REQUIRE(getContainerChunk(container, 0, 2) == std::vector<int>{0,1,2,3,4});
            REQUIRE(getContainerChunk(container, 1, 2) == std::vector<int>{5,6,7,8,9});
        }
        SECTION("Uneven Chunks") {
            REQUIRE(getContainerChunk(container, 0, 3) == std::vector<int>{0,1,2,3});
            REQUIRE(getContainerChunk(container, 1, 3) == std::vector<int>{4,5,6});
            REQUIRE(getContainerChunk(container, 2, 3) == std::vector<int>{7,8,9});
        }
        SECTION("Uneven Chunks 2") {
            REQUIRE(getContainerChunk(container, 0, 7) == std::vector<int>{0,1});
            REQUIRE(getContainerChunk(container, 1, 7) == std::vector<int>{2,3});
            REQUIRE(getContainerChunk(container, 2, 7) == std::vector<int>{4,5});
            REQUIRE(getContainerChunk(container, 3, 7) == std::vector<int>{6});
            REQUIRE(getContainerChunk(container, 4, 7) == std::vector<int>{7});
            REQUIRE(getContainerChunk(container, 5, 7) == std::vector<int>{8});
            REQUIRE(getContainerChunk(container, 6, 7) == std::vector<int>{9});
        }
        SECTION("Single Chunk") {
            REQUIRE(getContainerChunk(container, 0, 1) == std::vector<int>{0,1,2,3,4,5,6,7,8,9});
        }
    }
    SECTION("Unordered Map") {
        std::unordered_map<std::string, int> container;
        container.insert({"zero", 0});
        container.insert({"one", 1});
        container.insert({"two", 2});
        container.insert({"three", 3});
        container.insert({"four", 4});
        container.insert({"five", 5});
        container.insert({"six", 6});
        container.insert({"seven", 7});
        container.insert({"eight", 8});
        container.insert({"nine", 9});

        SECTION("More Threads Than Elements") {
            REQUIRE(getContainerChunk(container, 0, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 1, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 2, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 3, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 4, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 5, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 6, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 7, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 8, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 9, 11).size() == 1);
            REQUIRE(getContainerChunk(container, 10, 11).size() == 0);
            auto thrown = false;
            try {
                getContainerChunk(container, 11, 11);
            } catch (std::exception const & e) {
                thrown = true;
            }
            REQUIRE(thrown);
            // All chunks differnt
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 1, 11).begin()->first) == getContainerChunk(container, 0, 11).end());
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 2, 11).begin()->first) == getContainerChunk(container, 0, 11).end());
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 3, 11).begin()->first) == getContainerChunk(container, 0, 11).end());
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 4, 11).begin()->first) == getContainerChunk(container, 0, 11).end());
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 5, 11).begin()->first) == getContainerChunk(container, 0, 11).end());
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 6, 11).begin()->first) == getContainerChunk(container, 0, 11).end());
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 7, 11).begin()->first) == getContainerChunk(container, 0, 11).end());
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 8, 11).begin()->first) == getContainerChunk(container, 0, 11).end());
            REQUIRE(getContainerChunk(container, 0, 11).find(getContainerChunk(container, 9, 11).begin()->first) == getContainerChunk(container, 0, 11).end());

            REQUIRE(getContainerChunk(container, 1, 11).find(getContainerChunk(container, 2, 11).begin()->first) == getContainerChunk(container, 1, 11).end());
            REQUIRE(getContainerChunk(container, 1, 11).find(getContainerChunk(container, 3, 11).begin()->first) == getContainerChunk(container, 1, 11).end());
            REQUIRE(getContainerChunk(container, 1, 11).find(getContainerChunk(container, 4, 11).begin()->first) == getContainerChunk(container, 1, 11).end());
            REQUIRE(getContainerChunk(container, 1, 11).find(getContainerChunk(container, 5, 11).begin()->first) == getContainerChunk(container, 1, 11).end());
            REQUIRE(getContainerChunk(container, 1, 11).find(getContainerChunk(container, 6, 11).begin()->first) == getContainerChunk(container, 1, 11).end());
            REQUIRE(getContainerChunk(container, 1, 11).find(getContainerChunk(container, 7, 11).begin()->first) == getContainerChunk(container, 1, 11).end());
            REQUIRE(getContainerChunk(container, 1, 11).find(getContainerChunk(container, 8, 11).begin()->first) == getContainerChunk(container, 1, 11).end());
            REQUIRE(getContainerChunk(container, 1, 11).find(getContainerChunk(container, 9, 11).begin()->first) == getContainerChunk(container, 1, 11).end());

            REQUIRE(getContainerChunk(container, 3, 11).find(getContainerChunk(container, 4, 11).begin()->first) == getContainerChunk(container, 3, 11).end());
            REQUIRE(getContainerChunk(container, 3, 11).find(getContainerChunk(container, 5, 11).begin()->first) == getContainerChunk(container, 3, 11).end());
            REQUIRE(getContainerChunk(container, 3, 11).find(getContainerChunk(container, 6, 11).begin()->first) == getContainerChunk(container, 3, 11).end());
            REQUIRE(getContainerChunk(container, 3, 11).find(getContainerChunk(container, 7, 11).begin()->first) == getContainerChunk(container, 3, 11).end());
            REQUIRE(getContainerChunk(container, 3, 11).find(getContainerChunk(container, 8, 11).begin()->first) == getContainerChunk(container, 3, 11).end());
            REQUIRE(getContainerChunk(container, 3, 11).find(getContainerChunk(container, 9, 11).begin()->first) == getContainerChunk(container, 3, 11).end());

            REQUIRE(getContainerChunk(container, 4, 11).find(getContainerChunk(container, 5, 11).begin()->first) == getContainerChunk(container, 4, 11).end());
            REQUIRE(getContainerChunk(container, 4, 11).find(getContainerChunk(container, 6, 11).begin()->first) == getContainerChunk(container, 4, 11).end());
            REQUIRE(getContainerChunk(container, 4, 11).find(getContainerChunk(container, 7, 11).begin()->first) == getContainerChunk(container, 4, 11).end());
            REQUIRE(getContainerChunk(container, 4, 11).find(getContainerChunk(container, 8, 11).begin()->first) == getContainerChunk(container, 4, 11).end());
            REQUIRE(getContainerChunk(container, 4, 11).find(getContainerChunk(container, 9, 11).begin()->first) == getContainerChunk(container, 4, 11).end());

            REQUIRE(getContainerChunk(container, 5, 11).find(getContainerChunk(container, 6, 11).begin()->first) == getContainerChunk(container, 5, 11).end());
            REQUIRE(getContainerChunk(container, 5, 11).find(getContainerChunk(container, 7, 11).begin()->first) == getContainerChunk(container, 5, 11).end());
            REQUIRE(getContainerChunk(container, 5, 11).find(getContainerChunk(container, 8, 11).begin()->first) == getContainerChunk(container, 5, 11).end());
            REQUIRE(getContainerChunk(container, 5, 11).find(getContainerChunk(container, 9, 11).begin()->first) == getContainerChunk(container, 5, 11).end());

            REQUIRE(getContainerChunk(container, 6, 11).find(getContainerChunk(container, 7, 11).begin()->first) == getContainerChunk(container, 6, 11).end());
            REQUIRE(getContainerChunk(container, 6, 11).find(getContainerChunk(container, 8, 11).begin()->first) == getContainerChunk(container, 6, 11).end());
            REQUIRE(getContainerChunk(container, 6, 11).find(getContainerChunk(container, 9, 11).begin()->first) == getContainerChunk(container, 6, 11).end());

            REQUIRE(getContainerChunk(container, 7, 11).find(getContainerChunk(container, 8, 11).begin()->first) == getContainerChunk(container, 7, 11).end());
            REQUIRE(getContainerChunk(container, 7, 11).find(getContainerChunk(container, 9, 11).begin()->first) == getContainerChunk(container, 7, 11).end());

            REQUIRE(getContainerChunk(container, 8, 11).find(getContainerChunk(container, 9, 11).begin()->first) == getContainerChunk(container, 8, 11).end());
        }
        SECTION("Even Chunks") {
            REQUIRE(getContainerChunk(container, 0, 2).size() == 5);
            REQUIRE(getContainerChunk(container, 1, 2).size() == 5);
            auto chunk0 = getContainerChunk(container, 0, 2);
            auto chunk1 = getContainerChunk(container, 1, 2);
            for (auto&& elem : chunk0) {
                REQUIRE(chunk1.find(elem.first) == chunk1.end());
                REQUIRE(container.find(elem.first) != container.end());
            }
            for (auto&& elem : chunk1) {
                REQUIRE(chunk0.find(elem.first) == chunk0.end());
                REQUIRE(container.find(elem.first) != container.end());
            }

        }
        SECTION("Uneven Chunks") {
            REQUIRE(getContainerChunk(container, 0, 3).size() == 4);
            REQUIRE(getContainerChunk(container, 1, 3).size() == 3);
            REQUIRE(getContainerChunk(container, 2, 3).size() == 3);
            auto chunk0 = getContainerChunk(container, 0, 3);
            auto chunk1 = getContainerChunk(container, 1, 3);
            auto chunk2 = getContainerChunk(container, 2, 3);
            for (auto&& elem : chunk0) {
                REQUIRE(chunk1.find(elem.first) == chunk1.end());
                REQUIRE(chunk2.find(elem.first) == chunk2.end());
                REQUIRE(container.find(elem.first) != container.end());
            }
            for (auto&& elem : chunk1) {
                REQUIRE(chunk0.find(elem.first) == chunk0.end());
                REQUIRE(chunk2.find(elem.first) == chunk2.end());
                REQUIRE(container.find(elem.first) != container.end());
            }
            for (auto&& elem : chunk2) {
                REQUIRE(chunk0.find(elem.first) == chunk0.end());
                REQUIRE(chunk1.find(elem.first) == chunk1.end());
                REQUIRE(container.find(elem.first) != container.end());
            }
        }
        SECTION("Uneven Chunks 2") {
            REQUIRE(getContainerChunk(container, 0, 7).size() == 2);
            REQUIRE(getContainerChunk(container, 1, 7).size() == 2);
            REQUIRE(getContainerChunk(container, 2, 7).size() == 2);
            REQUIRE(getContainerChunk(container, 3, 7).size() == 1);
            REQUIRE(getContainerChunk(container, 4, 7).size() == 1);
            REQUIRE(getContainerChunk(container, 5, 7).size() == 1);
            REQUIRE(getContainerChunk(container, 6, 7).size() == 1);
        }
        SECTION("Single Chunk") {
            REQUIRE(getContainerChunk(container, 0, 1) == container);
        }
    }
}



TEST_CASE("Test containerChunkInsert") {
    std::cout << "[INFO] -- [TEST CASE] -- Test containerChunkInsert" << std::endl;
    auto array = std::array<int, 2>();
    auto map = std::map<int, int>();
    auto set = std::set<int>();
    auto unorderedMap = std::unordered_map<int, int>();
    auto unorderedSet = std::unordered_set<int>();
    auto vector = std::vector<int>();

    SECTION("Insert From Array") {
        array[0] = 0;
        array[1] = 1;
        SECTION("Array Batch") {
            containerChunkInsert<std::set<int>, std::array<int,2>>(set, array.begin(), array.end());
            REQUIRE(set == std::set<int>{0,1});

            containerChunkInsert<std::unordered_set<int>, std::array<int,2>>(unorderedSet, array.begin(), array.end());
            REQUIRE(unorderedSet.find(0) != unorderedSet.end());
            REQUIRE(unorderedSet.find(1) != unorderedSet.end());
            REQUIRE(unorderedSet.size() == 2);

            containerChunkInsert<std::vector<int>, std::array<int,2>>(vector, array.begin(), array.end());
            REQUIRE(vector == std::vector<int>{0,1});
        }
        SECTION("Array Single") {
            auto it = array.begin();
            containerChunkInsert<std::set<int>, std::array<int,2>>(set, it);
            containerChunkInsert<std::set<int>, std::array<int,2>>(set, ++it);
            REQUIRE(set == std::set<int>{0,1});

            it = array.begin();
            containerChunkInsert<std::unordered_set<int>, std::array<int,2>>(unorderedSet, it);
            containerChunkInsert<std::unordered_set<int>, std::array<int,2>>(unorderedSet, ++it);
            REQUIRE(unorderedSet.find(0) != unorderedSet.end());
            REQUIRE(unorderedSet.find(1) != unorderedSet.end());
            REQUIRE(unorderedSet.size() == 2);

            it = array.begin();
            containerChunkInsert<std::vector<int>, std::array<int,2>>(vector, it);
            containerChunkInsert<std::vector<int>, std::array<int,2>>(vector, ++it);
            REQUIRE(vector == std::vector<int>{0,1});
        }
    }

    SECTION("Insert From Map") {
        map.insert({0,10});
        map.insert({1,11});
        map.insert({2,12});
        std::map<int,int> target;
        std::unordered_map<int,int> unorderedTarget;
        SECTION("Map Batch") {
            containerChunkInsert<std::map<int,int>,std::map<int,int>>(target, map.begin(), map.end());
            REQUIRE(target == map);

            containerChunkInsert<std::unordered_map<int,int>,std::map<int,int>>(unorderedTarget, map.begin(), map.end());
            REQUIRE(unorderedTarget.at(0) == 10);
            REQUIRE(unorderedTarget.at(1) == 11);
            REQUIRE(unorderedTarget.at(2) == 12);
            REQUIRE(unorderedTarget.size() == 3);
        }
        SECTION("Map Single") {
            auto it = map.begin();
            containerChunkInsert<std::map<int,int>,std::map<int,int>>(target, it);
            containerChunkInsert<std::map<int,int>,std::map<int,int>>(target, ++it);
            containerChunkInsert<std::map<int,int>,std::map<int,int>>(target, ++it);
            REQUIRE(target == map);

            it = map.begin();
            containerChunkInsert<std::unordered_map<int,int>,std::map<int,int>>(unorderedTarget, it);
            containerChunkInsert<std::unordered_map<int,int>,std::map<int,int>>(unorderedTarget, ++it);
            containerChunkInsert<std::unordered_map<int,int>,std::map<int,int>>(unorderedTarget, ++it);
            REQUIRE(unorderedTarget.at(0) == 10);
            REQUIRE(unorderedTarget.at(1) == 11);
            REQUIRE(unorderedTarget.at(2) == 12);
            REQUIRE(unorderedTarget.size() == 3);
        }

        // insert from unoredered_map
        unorderedMap.insert({0,10});
        unorderedMap.insert({1,11});
        unorderedMap.insert({2,12});
        SECTION("UMap Batch") {
            containerChunkInsert<std::map<int,int>,std::unordered_map<int,int>>(target, unorderedMap.begin(), unorderedMap.end());
            REQUIRE(target == map);

            containerChunkInsert<std::unordered_map<int,int>,std::unordered_map<int,int>>(unorderedTarget, unorderedMap.begin(), unorderedMap.end());
            REQUIRE(unorderedTarget.at(0) == 10);
            REQUIRE(unorderedTarget.at(1) == 11);
            REQUIRE(unorderedTarget.at(2) == 12);
            REQUIRE(unorderedTarget.size() == 3);
        }
        SECTION("UMap Single") {
            auto it = unorderedMap.begin();
            containerChunkInsert<std::map<int,int>,std::unordered_map<int,int>>(target, it);
            containerChunkInsert<std::map<int,int>,std::unordered_map<int,int>>(target, ++it);
            containerChunkInsert<std::map<int,int>,std::unordered_map<int,int>>(target, ++it);
            REQUIRE(target == map);

            it = unorderedMap.begin();
            containerChunkInsert<std::unordered_map<int,int>,std::unordered_map<int,int>>(unorderedTarget, it);
            containerChunkInsert<std::unordered_map<int,int>,std::unordered_map<int,int>>(unorderedTarget, ++it);
            containerChunkInsert<std::unordered_map<int,int>,std::unordered_map<int,int>>(unorderedTarget, ++it);
            REQUIRE(unorderedTarget.at(0) == 10);
            REQUIRE(unorderedTarget.at(1) == 11);
            REQUIRE(unorderedTarget.at(2) == 12);
            REQUIRE(unorderedTarget.size() == 3);
        }
    }

    SECTION("Insert From Set") {
        set.emplace(0);
        set.emplace(1);
        set.emplace(2);
        std::set<int> targetSet;
        std::unordered_set<int> targetUnorderedSet;
        std::vector<int> targetVector;
        SECTION("Set Batch") {
            containerChunkInsert<std::set<int>, std::set<int>>(targetSet, set.begin(), set.end());
            REQUIRE(targetSet == set);

            containerChunkInsert<std::unordered_set<int>, std::set<int>>(targetUnorderedSet, set.begin(), set.end());
            REQUIRE(targetUnorderedSet.find(0) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(1) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(2) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.size() == 3);

            containerChunkInsert<std::vector<int>, std::set<int>>(targetVector, set.begin(), set.end());
            REQUIRE(targetVector == std::vector<int>{0,1,2});
        }
        SECTION("Set Single") {
            auto it = set.begin();
            containerChunkInsert<std::set<int>, std::set<int>>(targetSet, it);
            containerChunkInsert<std::set<int>, std::set<int>>(targetSet, ++it);
            containerChunkInsert<std::set<int>, std::set<int>>(targetSet, ++it);
            REQUIRE(targetSet == set);

            it = set.begin();
            containerChunkInsert<std::unordered_set<int>, std::set<int>>(targetUnorderedSet, it);
            containerChunkInsert<std::unordered_set<int>, std::set<int>>(targetUnorderedSet, ++it);
            containerChunkInsert<std::unordered_set<int>, std::set<int>>(targetUnorderedSet, ++it);
            REQUIRE(targetUnorderedSet.find(0) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(1) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(2) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.size() == 3);

            it = set.begin();
            containerChunkInsert<std::vector<int>, std::set<int>>(targetVector, it);
            containerChunkInsert<std::vector<int>, std::set<int>>(targetVector, ++it);
            containerChunkInsert<std::vector<int>, std::set<int>>(targetVector, ++it);
            REQUIRE(targetVector == std::vector<int>{0,1,2});
        }

        // insert from unordered_set
        unorderedSet.emplace(0);
        unorderedSet.emplace(1);
        unorderedSet.emplace(2);
        SECTION("USet Batch") {
            containerChunkInsert<std::set<int>, std::unordered_set<int>>(targetSet, unorderedSet.begin(), unorderedSet.end());
            REQUIRE(targetSet == set);

            containerChunkInsert<std::unordered_set<int>, std::unordered_set<int>>(targetUnorderedSet, unorderedSet.begin(), unorderedSet.end());
            REQUIRE(targetUnorderedSet.find(0) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(1) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(2) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.size() == 3);

            containerChunkInsert<std::vector<int>, std::unordered_set<int>>(targetVector, unorderedSet.begin(), unorderedSet.end());
            std::sort(targetVector.begin(), targetVector.end());
            REQUIRE(targetVector == std::vector<int>{0,1,2});
        }
        SECTION("USet Single") {
            auto it = unorderedSet.begin();
            containerChunkInsert<std::set<int>, std::unordered_set<int>>(targetSet, it);
            containerChunkInsert<std::set<int>, std::unordered_set<int>>(targetSet, ++it);
            containerChunkInsert<std::set<int>, std::unordered_set<int>>(targetSet, ++it);
            REQUIRE(targetSet == set);

            it = unorderedSet.begin();
            containerChunkInsert<std::unordered_set<int>, std::unordered_set<int>>(targetUnorderedSet, it);
            containerChunkInsert<std::unordered_set<int>, std::unordered_set<int>>(targetUnorderedSet, ++it);
            containerChunkInsert<std::unordered_set<int>, std::unordered_set<int>>(targetUnorderedSet, ++it);
            REQUIRE(targetUnorderedSet.find(0) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(1) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(2) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.size() == 3);

            it = unorderedSet.begin();
            containerChunkInsert<std::vector<int>, std::unordered_set<int>>(targetVector, it);
            containerChunkInsert<std::vector<int>, std::unordered_set<int>>(targetVector, ++it);
            containerChunkInsert<std::vector<int>, std::unordered_set<int>>(targetVector, ++it);
            std::sort(targetVector.begin(), targetVector.end());
            REQUIRE(targetVector == std::vector<int>{0,1,2});
        }
    }

    SECTION("Insert From Vector") {
        vector.emplace_back(0);
        vector.emplace_back(1);
        vector.emplace_back(2);
        std::set<int> targetSet;
        std::unordered_set<int> targetUnorderedSet;
        std::vector<int> targetVector;
        SECTION("Vector Batch") {
            containerChunkInsert<std::set<int>, std::vector<int>>(targetSet, vector.begin(), vector.end());
            REQUIRE(targetSet == std::set<int>{0,1,2});

            containerChunkInsert<std::unordered_set<int>, std::vector<int>>(targetUnorderedSet, vector.begin(), vector.end());
            REQUIRE(targetUnorderedSet.find(0) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(1) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(2) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.size() == 3);

            containerChunkInsert<std::vector<int>, std::vector<int>>(targetVector, vector.begin(), vector.end());
            REQUIRE(targetVector == vector);
        }
        SECTION("Vector Single") {
            auto it = vector.begin();
            containerChunkInsert<std::set<int>, std::vector<int>>(targetSet, it);
            containerChunkInsert<std::set<int>, std::vector<int>>(targetSet, ++it);
            containerChunkInsert<std::set<int>, std::vector<int>>(targetSet, ++it);
            REQUIRE(targetSet == std::set<int>{0,1,2});

            it = vector.begin();
            containerChunkInsert<std::unordered_set<int>, std::vector<int>>(targetUnorderedSet, it);
            containerChunkInsert<std::unordered_set<int>, std::vector<int>>(targetUnorderedSet, ++it);
            containerChunkInsert<std::unordered_set<int>, std::vector<int>>(targetUnorderedSet, ++it);
            REQUIRE(targetUnorderedSet.find(0) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(1) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.find(2) != targetUnorderedSet.end());
            REQUIRE(targetUnorderedSet.size() == 3);

            it = vector.begin();
            containerChunkInsert<std::vector<int>, std::vector<int>>(targetVector, it);
            containerChunkInsert<std::vector<int>, std::vector<int>>(targetVector, ++it);
            containerChunkInsert<std::vector<int>, std::vector<int>>(targetVector, ++it);
            REQUIRE(targetVector == vector);
        }
    }
}
