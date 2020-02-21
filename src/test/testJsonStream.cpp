#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "catch2/catch.hpp"
#include "../JsonStream.h"

TEST_CASE("JsonStream") {
    SECTION("Test JsonValue") {
        auto str = std::string("\"tr\"y this\"");
        auto val = JsonValue(str);
        REQUIRE(val.value() == "\"tr\\\"y this\"");

        auto doubleNum = JsonValue(1.);
        REQUIRE(doubleNum.value() == std::to_string(1.));
        auto intNum = JsonValue(1);
        REQUIRE(intNum.value() == "1");

        std::vector<int> innerVector{1,2,3,4,5};
        std::vector<std::vector<int>> outerVector;
        outerVector.emplace_back(innerVector);
        outerVector.emplace_back(innerVector);
        auto jvector = JsonValue(outerVector);
        REQUIRE(jvector.value() == "[[1,2,3,4,5],[1,2,3,4,5]]");
        REQUIRE(JsonValue(std::vector<int>{}).value() == "[]");

        auto set = std::set<int>({1,2,3,4,5});
        auto jset = JsonValue(set);
        REQUIRE(jset.value() == "[1,2,3,4,5]");
        REQUIRE(JsonValue(std::set<int>{}).value() == "[]");

        auto map = std::map<int,std::string>({{1, "1"}, {2,"2"}});
        auto jmap = JsonValue(map);
        REQUIRE(jmap.value() == "{\"1\":\"1\",\"2\":\"2\"}");
        REQUIRE(JsonValue(std::map<int,int>{}).value() == "{}");

        REQUIRE(JsonValue(false).value() == "false");
        REQUIRE(JsonValue(true).value() == "true");

        REQUIRE(JsonValue("something", JsonValue::stringIsValidJson).value() == "something");
    }

    SECTION("Test Streaming") {
        auto outstream = std::stringstream(std::ios_base::out);
        SECTION("JsonStreamArray") {
            auto jstream = JsonStreamArray(outstream);
            auto array = std::array<int, 2>{4,2};
            auto jsonArray = JsonValue(array);
            jstream << array << array;
            jstream.addValue(array);
            jstream << jsonArray;
            jstream.addValue(jsonArray);
            jstream.close();
            auto string = outstream.str();
            REQUIRE(string == "[[4,2],[4,2],[4,2],[4,2],[4,2]]");
        }
        SECTION("JsonStreamDict") {
            auto jstream = JsonStreamDict(outstream);
            std::map<int, int> map;
            map[1] = 10;
            map[2] = 20;
            map[3] = 30;
            for (auto&& elem : map) {
                jstream.addValue(std::to_string(elem.first), elem.second);
            }
            jstream.addValue("somethingCompletelyDifferent", true);
            JsonValue innerDict{std::map<int,bool>{{0, false}}};
            jstream.addValue("innerDict", innerDict);
            jstream.close();
            auto string = outstream.str();
            REQUIRE(string == "{\"1\":10,\"2\":20,\"3\":30,\"somethingCompletelyDifferent\":true,\"innerDict\":{\"0\":false}}");
        }
    }
}
