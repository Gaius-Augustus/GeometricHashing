#ifndef GEOMETRICHASHING_H
#define GEOMETRICHASHING_H

#include <memory>
#include <ostream>
#include <sstream>

#include "json/json.hpp"
#include "Configuration.h"
#include "Cubeset.h"
#include "IdentifierMapping.h"
#include "JsonStream.h"
#include "Linkset.h"
#include "SeedMapGeneral.h"
#include "Timestep.h"

template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
inline void geometricHashing(std::shared_ptr<IdentifierMapping const> idMap,
                             std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> seedMap,
                             std::shared_ptr<FastaCollection const> fastaCollection,
                             std::shared_ptr<Configuration const> config,
                             bool cubeOutput = true) {
    seedMap->clearMatches();
    auto linkset = std::make_shared<Linkset>(idMap,
                                             config->matchLimit(),
                                             config->matchLimitDiscardSeeds(),
                                             config->occurrencePerGenomeMax(),
                                             config->occurrencePerGenomeMin());
    Timestep tsLinkset("Creating Linkset");
    std::bitset<16> fullyDimensional;
    for (auto&& elem : seedMap->seedMap()) {
        for (auto&& occurrences : elem.second) {
            for (auto&& occurrence : occurrences) {
                fullyDimensional.set(occurrence.genome());
            }
            // only create fully-dimensional links
            if (fullyDimensional.count() == idMap->numGenomes()) {
                linkset->createLinks(occurrences);
            }
            fullyDimensional.reset();
        }
    }
    std::cout << "Created " << linkset->numLinks() << " Links" << std::endl;
    tsLinkset.endAndPrint(Timestep::minutes);
    // Attention: Links now use _exact_ positions instead of tile IDs
    Timestep tsCubeset("Creating Cubeset");
    Cubeset cubeset(linkset, fastaCollection, config);
    std::cout << "Created " << cubeset.cubeMap().size() << " Cubes" << std::endl;
    tsCubeset.endAndPrint(Timestep::seconds);

    Timestep tsMatches("Extracting Matches from Cubes");
    seedMap->clearMatches();

    if (cubeOutput) { outputCubes(cubeset, seedMap, config); }

    for (auto&& elem : cubeset.scoreToCube()) {
        auto score = elem.first;
        auto& cubes = elem.second;
        if (score >= static_cast<double>(config->cubeScoreThreshold())) {
            for (auto&& cubeptr : cubes) {
                auto& links = cubeset.cubeMap().at(cubeptr);
                for (auto&& link : links) {
                    seedMap->createMatchesFromOccurrences(link->occurrence());
                }
            }
        }
    }
    std::cout << "Created " << seedMap->numMatches() << " matches" << std::endl;
    tsMatches.endAndPrint(Timestep::seconds);
}



template<typename TwoBitKmerDataType, typename TwoBitSeedDataType>
inline void outputCubes(Cubeset const & cubeset,
                        std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> seedMap,
                        std::shared_ptr<Configuration const> config) {
    auto idMap = seedMap->idMap();
    auto outpath = config->output();
    outpath.replace_extension(".cubes.json");
    auto os = std::ofstream(outpath);
    auto jstream = JsonStreamDict(os, true);
    for (auto&& elem : cubeset.scoreToCube()) {
        auto score = elem.first;
        auto& cubes = elem.second;
        if (score >= config->cubeScoreThreshold()) {
            for (auto&& cubeptr : cubes) {
                // create custom key to identify a cube
                std::vector<JsonValue> keyVector;
                for (auto&& td : cubeptr->tiledistance()) {
                    auto keyParts = std::array<std::string, 4>{idMap->queryGenomeName(td.genome()),
                                                               idMap->querySequenceName(td.sequence()),
                                                               std::to_string(td.distance()),
                                                               std::to_string(td.reverse())};
                    keyVector.emplace_back(keyParts);
                }
                auto cubeKey = JsonValue(keyVector);
                // collect links in cube and their counts
                std::map<std::string, double> linkToCount;
                auto& links = cubeset.cubeMap().at(cubeptr);
                for (auto&& link : links) {
                    std::vector<JsonValue> occurrences;
                    for (auto&& occ : link->occurrence()) {
                        occurrences.emplace_back(occ.toJsonValue(seedMap->span(),
                                                 seedMap->idMap()));
                    }
                    linkToCount.emplace(JsonValue(occurrences).value(),
                                        cubeset.underlyingLinkset()->linkCount(link));
                }
                linkToCount.emplace("score", score);    // also include cube score
                linkToCount.emplace("rawCount", links.size()); // and raw link count
                jstream.addValue(cubeKey.value(), linkToCount);
            }
        }
    }
    jstream.close();
}

#endif // GEOMETRICHASHING_H
