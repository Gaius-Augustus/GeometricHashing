/* External Libraries
 * ~~~~~~~~~~~~~~~~~~
 * - Boost 1.70.0
 * - cxx-prettyprint (https://github.com/louisdx/cxx-prettyprint) 
 * - Catch 2 (https://github.com/catchorg/Catch2)
 * - JSON for Modern C++ version 3.1.2 (https://github.com/nlohmann/json)
 * ~~~~~~~~~~~~~~~~~~ */

/* Standard include ordering: Related header, (blank line), C library, C++ library, (blank line), other libraries' .h, your project's .h.
 * Inside the blocks: alphabetical order */

#include <climits>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>

#include "Configuration.h"
#include "Cubeset.h"
#include "ExtractSeeds.h"
#include "FastaCollection.h"
#include "GeometricHashing.h"
#include "JsonStream.h"
#include "Linkset.h"
#include "MemoryMonitor.h"
#include "SeedMapContiguous.h"
#include "SeedMapSpaced.h"
#include "Timestep.h"
#include "TwoBitKmer.h"

namespace fs = std::experimental::filesystem;
namespace po = boost::program_options;



template <typename TwoBitKmerDataType, typename TwoBitSeedDataType,
          template<typename, typename> typename SeedMapChild>
auto runSeedExtractionAndFiltering(std::shared_ptr<Configuration const> config) {
    Timestep tsRun{"Run"};
    auto jsonstream = std::stringstream(std::ios_base::out);    // stream json of run info into string
    JsonStreamDict runInfo(jsonstream);

    // prepare
    std::ofstream os(config->output());
    if (!os.good()) { throw std::runtime_error("[ERROR] -- Cannot write to outfile"); }

    auto fastaCollection = std::make_shared<FastaCollection>(config);   // create all FastaRepresentations, including the artificial ones
    std::cout << "[INFO] -- Total number of sequences: " << fastaCollection->numSequences() << std::endl << std::endl;
    if (config->outputArtificialSequences().string() != "") {
        std::ofstream os(config->outputArtificialSequences());
        for (auto&& rep : fastaCollection->collection()) {
            rep.second.writeArtificialSequences(os);
        }
    }

    auto idMap = std::make_shared<IdentifierMapping>(config->genome1());
    idMap->queryGenomeID(config->genome2());
    fastaCollection->populateIdentifierMappingFromFastaCollection(*idMap);
    std::shared_ptr<SpacedSeedMaskCollection const> masks;
    if (config->masks().size() > 0) {
        masks = std::make_shared<SpacedSeedMaskCollection const>(config->masks());
    } else if (config->optimalSeed()) {
        masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                 SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
    } else {
        masks = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(config->weight()),
                                                                 SpacedSeedMaskCollection::Span(config->span()),
                                                                 SpacedSeedMaskCollection::SeedSetSize(config->seedSetSize()));
    }

    std::cout << "[INFO] -- Masks used: " << *masks << std::endl << std::endl;
    runInfo.addValue("masks", masks->masksAsString());

    // run seed extraction and initial match creation
    auto seedMap = std::make_shared<SeedMapChild<TwoBitKmerDataType,
                                                 TwoBitSeedDataType>>(config, idMap, masks);
    Timestep tsDirectExtract("Extracting and filtering seeds from input data");
    // last parameter: skipMatchCreation, if --geometric-hashing is set, no need for matches (experimental)
    auto extract = ExtractSeeds<TwoBitKmerDataType,
                                TwoBitSeedDataType>(fastaCollection,
                                                    seedMap,
                                                    config->nThreads(), config->quiet(), config->performGeometricHashing());
    tsDirectExtract.endAndPrint(Timestep::seconds);

    // print some statistics stuff
    std::map<size_t, std::map<size_t, size_t>> directCounts;
    for (auto&& match : seedMap->matches()) {
        auto i = match.first().genome();
        auto j = match.second().genome();
        ++(directCounts[i][j]);
    }
    for (auto&& elem1 : directCounts) {
        for (auto&& elem2 : elem1.second) {
            std::cout << "[STATISTICSTAG001] direct matches " << idMap->queryGenomeName(elem1.first) << "-" << idMap->queryGenomeName(elem2.first) << " " << elem2.second << std::endl;
            std::string key{"directMatches_"};
            key += idMap->queryGenomeName(elem1.first) + "-" + idMap->queryGenomeName(elem2.first);
            runInfo.addValue(key, elem2.second);
        }
    }

    // run advanced filtering
    if (config->performGeometricHashing()) {
        // this includes M4 filtering inside Cubes, M5 filter makes no sense here
        std::cout << "[INFO] -- Running with GeometricHashing Filter" << std::endl;
        geometricHashing<TwoBitKmerDataType, TwoBitSeedDataType>(idMap, seedMap, config);
    } else {
        // M4 and M5 may both be applied (first M5 to get more matches, then M4 to filter them)

        // save memory
        seedMap->clearSeedMap();
        if (config->performDiagonalFiltering()) {
            std::cout << "[INFO] -- Running with Diagonal Matches Filter" << std::endl;
            seedMap->applyDiagonalMatchesFilter(std::make_shared<std::array<size_t, 4>>());
        }
    }

    // save memory
    seedMap->clearSeedMap();

    std::cout << "Writing matches to output file..." << std::endl;
    seedMap->output(os);
    os.close();

    tsRun.end();
    runInfo.addValue("runtime", tsRun.elapsed(Timestep::minutes));
    runInfo.close();
    return JsonValue(jsonstream.str(), JsonValue::stringIsValidJson);
}



//! Main
/*! Intermediate outputs and chrono stuff is for debugging purposes and not crucial
 * for this program. */
int main(int argc, char * argv[]) {
    // parse cmdline options
    auto config = std::make_shared<Configuration const>(argc, argv);

    // Print info on parameters
    std::cout << std::endl << "Run Seed Finding" << std::endl
              <<              "~~~~~~~~~~~~~~~~" << std::endl << std::endl;
    std::cout << "Processing files" << std::endl;
    for (auto&& file : config->inputFiles()) { std::cout << "\t" << file << std::endl; }
    std::cout << "With " << *config << std::endl;
    std::cout << std::endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Timestep tsMain("Running Program");

    auto mm = MM::MemoryMonitor();
    std::cout << mm << std::endl << std::endl;
    std::cout << "Starting Program" << std::endl;

    JsonValue runInfo;
    if (config->weight() == config->span() && !config->performDiagonalFiltering()) {
        std::cout << "[INFO] -- Running with Contiguous Seeds" << std::endl;
        if (config->span() <= 29) {
            runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataShort, TwoBitKmerDataShort, SeedMapContiguous>(config);
        } else if (config->span() <= 61) {
            runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataMedium, TwoBitKmerDataMedium, SeedMapContiguous>(config);
        } else {
            runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataLong, TwoBitKmerDataLong, SeedMapContiguous>(config);
        }
    } else {
        std::cout << "[INFO] -- Running with Spaced Seeds" << std::endl;
        if (config->span() <= 29) {
            runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataShort, TwoBitKmerDataShort, SeedMapSpaced>(config);
        } else if (config->span() <= 61) {
            if (config->weight() <= 29) {
                runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataMedium, TwoBitKmerDataShort, SeedMapSpaced>(config);
            } else {
                runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataMedium, TwoBitKmerDataMedium, SeedMapSpaced>(config);
            }
        } else {
            if (config->weight() <= 29) {
                runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataLong, TwoBitKmerDataShort, SeedMapSpaced>(config);
            } else if (config->weight() <= 61) {
                runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataLong, TwoBitKmerDataMedium, SeedMapSpaced>(config);
            } else {
                runInfo = runSeedExtractionAndFiltering<TwoBitKmerDataLong, TwoBitKmerDataLong, SeedMapSpaced>(config);
            }
        }
    }

    if (config->outputRunInformation().string() != "") {
        auto configJson = config->configJson();
        std::ofstream os(config->outputRunInformation());
        if (os.good()) {
            JsonStreamDict json{os};
            json.addValue("configuration", configJson);
            json.addValue("run", runInfo);
            json.close();
        } else {
            std::cerr << "[WARNING] -- Cannot write to run information output file '" << config->outputRunInformation() << "'";
        }
    }

    std::cout << "Finished Program" << std::endl;
    std::cout << mm << std::endl;
    tsMain.endAndPrint(Timestep::minutes);

    return 0;
}
