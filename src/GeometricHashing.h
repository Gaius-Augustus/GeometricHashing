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
inline void geometricHashing(std::shared_ptr<IdentifierMapping> idMap,
                             std::shared_ptr<SeedMapGeneral<TwoBitKmerDataType, TwoBitSeedDataType>> seedMap,
                             std::shared_ptr<Configuration const> config) {
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
    Cubeset cubeset(linkset, config);
    std::cout << "Created " << cubeset.cubeMap().size() << " Cubes" << std::endl;
    tsCubeset.endAndPrint(Timestep::seconds);

    Timestep tsMatches("Extracting Matches from Cubes");
    seedMap->clearMatches();

    outputCubes(cubeset, seedMap, config);

    for (auto&& elem : cubeset.scoreToCube()) {
        auto score = elem.first;
        auto& cubes = elem.second;
        if (score >= config->cubeScoreThreshold()) {
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
    auto jstream = JsonStreamDict(os);
    auto first = true;
    for (auto&& elem : cubeset.scoreToCube()) {
        auto score = elem.first;
        auto& cubes = elem.second;
        if (score >= config->cubeScoreThreshold()) {
            for (auto&& cubeptr : cubes) {
                // create custom key to identify a cube
                std::ostringstream cubeKey;
                size_t i = 0;
                for (auto&& td : cubeptr->tiledistance()) {
                    cubeKey << "(" << idMap->queryGenomeName(td.genome()) << ", " << idMap->querySequenceName(td.sequence()) << ", " << std::to_string(td.distance()) << ", " << td.reverse() << ")";
                    ++i; if (i < cubeptr->tiledistance().size()) { cubeKey << ", "; }
                }
                // simluate an output json array and write the first two occurrences of each respective link to it
                std::ostringstream linkstr;
                auto linkstream = JsonStreamArray(linkstr);
                auto& links = cubeset.cubeMap().at(cubeptr);
                for (auto&& link : links) {
                    seedMap->appendMatchToOutput(linkstream, link->occurrence(0), link->occurrence(1));
                }
                linkstream.close();
                // add json array string as value to cube key
                if (!first) { jstream.ostream() << ",";} else { first = false; }
                jstream.ostream() << JsonValue(cubeKey.str()).value() << ":" << linkstr.str();
            }
        }
    }
    jstream.close();
}

#endif // GEOMETRICHASHING_H
