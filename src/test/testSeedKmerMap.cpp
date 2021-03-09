#include "catch2/catch.hpp"
#include "LoadTestdata.h"
#include "../SeedKmerMap.h"

TEST_CASE("Load Testdata into SeedKmerMap") {
    std::cout << "[INFO] -- [TEST CASE] -- Load Testdata into SeedKmerMap" << std::endl;

    auto valueEqual = [](std::vector<std::vector<KmerOccurrence>> lhs,
                         std::vector<std::vector<KmerOccurrence>> rhs) {
        if (lhs.size() == rhs.size()) {
            for (size_t i = 0; i < lhs.size(); ++i) {
                if (lhs[i].size() > 1) {
                    std::sort(lhs[i].begin(), lhs[i].end());
                    std::sort(rhs[i].begin(), rhs[i].end());
                }
                if (lhs[i] != rhs[i]) {
                    std::cout << "[INFO] -- value equal -- " << lhs << " != " << rhs << std::endl;
                    return false;
                }
            }
            return true;
        }
        std::cout << "[INFO] -- value equal -- " << lhs << " != " << rhs << std::endl;
        return false;
    };

    auto data = LoadTestdataMetagraph();
    std::cout << "[INFO] -- Executing with " << data.configGraph()->nThreads() << " threads" << std::endl;
    auto expectedSeedsGraph = data.expectedSeedsGraph(*data.idMapGraph());
    auto expectedLinksGraph = data.expectedLinksGraph(*data.idMapGraph());
    auto idMap = std::make_shared<IdentifierMapping>(data.configGraph()->genome1());
    *idMap = *data.idMapGraph(); // copy IDs

    SeedKmerMap<TwoBitKmerDataShort> map{data.configGraph(), idMap};
    map.extractSeeds();
    REQUIRE(map.size() == expectedSeedsGraph.size());
    for (auto&& elem : expectedSeedsGraph) {
        REQUIRE(map.map().find(elem.first) != map.map().end());
    }
    // direct map.begin() iteration takes very long
    for (auto&& elem : map.map()) {
        REQUIRE(expectedSeedsGraph.find(elem.first) != expectedSeedsGraph.end());
    }
    auto i = 0;
    for (auto&& elem : expectedSeedsGraph) {
        REQUIRE(valueEqual(map.at(elem.first), elem.second));
        ++i;
        if (i >= 1000) { break; } // takes very long
    }
    // parallel iteration
    auto callback = [&map,
                     &expectedSeedsGraph,
                     &valueEqual](typename SeedKmerMap<TwoBitKmerDataShort>::SeedKmerMapType::const_iterator it,
                                  typename SeedKmerMap<TwoBitKmerDataShort>::SeedKmerMapType::const_iterator end) {
        auto i = 0;
        for (; it != end; ++it) {
            REQUIRE(expectedSeedsGraph.find(it->first) != expectedSeedsGraph.end());
            REQUIRE(valueEqual(expectedSeedsGraph.at(it->first), map.at(it->first)));
            if (i >= 1000) { break; }
        }
    };
    map.parallelIteration(callback);
    REQUIRE(map.size() == expectedSeedsGraph.size());
    // cleanup and create links
    map.cleanup();
    Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan> observedLinkset{data.configGraph(), idMap};
    for (auto&& elem : map.map()) {
        auto occs = map.at(elem.first);
        for (size_t maskInd = 0; maskInd < data.configGraph()->seedSetSize(); ++maskInd) {
            auto& occv = occs.at(maskInd);
        //for (auto&& occv : occs) {
            observedLinkset.createLinks(occv, data.configGraph()->maskCollection()->span(maskInd));
        }
    }
    REQUIRE(hopscotchEqual(observedLinkset.linkset(), expectedLinksGraph));
}
