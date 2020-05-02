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

        str = "[[\"hg38\",\"hg38|chr18|52340197|1195707|+|80373285|1195707|ENSG00000187323|ENST00000442544\",\"0\",\"0\"],[\"mm10\",\"mm10|chr18|71253634|1097436|-|90702639|1097436|ENSMUSG00000060534\",\"-8\",\"0\"]]";
        val = JsonValue(str);
        REQUIRE(val.value() == "\"[[\\\"hg38\\\",\\\"hg38|chr18|52340197|1195707|+|80373285|1195707|ENSG00000187323|ENST00000442544\\\",\\\"0\\\",\\\"0\\\"],[\\\"mm10\\\",\\\"mm10|chr18|71253634|1097436|-|90702639|1097436|ENSMUSG00000060534\\\",\\\"-8\\\",\\\"0\\\"]]\"");

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

        // test some use cases to make sure output is always correctly formatted
        auto innerStr0 = std::string("somethingSomething");
        auto innerStrJson0 = JsonValue(innerStr0);
        REQUIRE(innerStrJson0 == JsonValue(innerStrJson0));
        auto innerStr1 = std::string("somethingElse");
        auto innerStrJson1 = JsonValue(innerStr1);
        auto innerStr2 = std::string("somethingDifferent");
        auto innerStrJson2 = JsonValue(innerStr2);
        auto innerStr3 = std::string("somethingCompletelyDifferent");
        auto innerStrJson3 = JsonValue(innerStr3);
        auto innerArray0 = std::array<std::string, 2>{innerStr0, innerStr1};
        auto innerJsonArray0 = std::array<JsonValue, 2>{innerStrJson0, innerStrJson1};
        auto innerArray1 = std::array<std::string, 2>{innerStr2, innerStr3};
        auto innerJsonArray1 = std::array<JsonValue, 2>{innerStrJson2, innerStrJson3};
        REQUIRE(JsonValue(innerArray0) == JsonValue(innerJsonArray0));
        REQUIRE(JsonValue(innerArray1) == JsonValue(innerJsonArray1));
        auto outerArray = std::array<std::array<std::string, 2>, 2>{innerArray0, innerArray1};
        auto outerJsonArray = std::array<JsonValue, 2>{JsonValue(innerJsonArray0), JsonValue(innerJsonArray1)};
        REQUIRE(JsonValue(outerArray) == JsonValue(outerJsonArray));
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
        SECTION("JsonStreamArray with newline") {
            auto jstream = JsonStreamArray(outstream, true);
            auto array = std::array<int, 2>{4,2};
            auto jsonArray = JsonValue(array);
            jstream << array << array;
            jstream.addValue(array);
            jstream << jsonArray;
            jstream.addValue(jsonArray);
            jstream.close();
            auto string = outstream.str();
            REQUIRE(string == "[[4,2],\n[4,2],\n[4,2],\n[4,2],\n[4,2]]");
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
        SECTION("JsonStreamDict with newline") {
            auto jstream = JsonStreamDict(outstream, true);
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
            REQUIRE(string == "{\"1\":10,\n\"2\":20,\n\"3\":30,\n\"somethingCompletelyDifferent\":true,\n\"innerDict\":{\"0\":false}}");
        }
    }
}
