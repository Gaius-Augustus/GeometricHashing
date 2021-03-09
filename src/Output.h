#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <memory>

#include "mabl3/JsonStream.h"
#include "Configuration.h"
#include "Cubeset.h"
#include "KmerOccurrence.h"
#include "regionTupleExtraction.h"

using namespace mabl3;

class Output {
public:
    //! c'tor
    /*! \param config Configuration object */
    Output(std::shared_ptr<Configuration const> config)
        : config_{config},
          cubeDict_{},
          cubeOs_{},
          linkOs_{config_->output()},
          linkArray_{linkOs_},
          mutexOutput_{},
          relevantCubeDict_{},
          relevantCubeOs_{},
          runInfo_{},
          runOs_{} {
        if (!linkOs_.good()) { throw std::runtime_error("[ERROR] -- Output -- Cannot write to outfile"); }
        // Write config to run information if demanded
        if (config->outputRunInformation().string() != "") {
            runOs_ = std::make_unique<std::ofstream>(config->outputRunInformation());
            if (runOs_->good()) {
                runInfo_ = std::make_unique<JsonStreamDict>(*runOs_, true);
                runInfo_->addValue("configuration", config->configJson());
                runInfo_->addValue("masks", config_->maskCollection()->masksAsString());
            } else {
                std::cerr << "[WARNING] -- Cannot write to run information output file '" << config->outputRunInformation() << "'" << std::endl;
            }
        }
        // Open cube stream if demanded
        if (config_->performGeometricHashing()) {
            auto outpath = config_->output();
            outpath.replace_extension(".cubes.json");
            cubeOs_ = std::make_unique<std::ofstream>(outpath);
            if (cubeOs_->good()) {
                cubeDict_ = std::make_unique<JsonStreamDict>(*cubeOs_, true);
            } else {
                std::cerr << "[WARNING] -- Cannot write cubes to output file '" << outpath << "'" << std::endl;
            }
        }
        // Open relevant cube stream if appropriate
        if (config_->preMaskCollection()) {
            auto outpath = config_->output();
            outpath.replace_extension(".relevantCubes.json");
            relevantCubeOs_ = std::make_unique<std::ofstream>(outpath);
            if (relevantCubeOs_->good()) {
                relevantCubeDict_ = std::make_unique<JsonStreamDict>(*relevantCubeOs_, true);
            } else {
                std::cerr << "[WARNING] -- Cannot write relevant cubes to output file '" << outpath << "'" << std::endl;
            }
        }
    }
    ~Output() {
        // close json streams before underlying ofstreams get closed
        if (cubeDict_) { cubeDict_->close(); }
        linkArray_.close();
        if (relevantCubeDict_) { relevantCubeDict_->close(); }
        if (runInfo_) { runInfo_->close(); }
    }

    //! Add a kay-value pair to run info
    template<typename T>
    void addRunInfo(std::string const & key, T const & value) {
        if (runInfo_) {
            std::unique_lock<std::mutex> outputLock(mutexOutput_);
            runInfo_->addValue(key, value);
            outputLock.unlock();
        }
    }
    //! Append a single Link to link output
    template <typename LinkType>
    void appendLinkToLinkOutput(LinkType const & link,
                                IdentifierMapping const & idMap) {
        std::unique_lock<std::mutex> outputLock(mutexOutput_);
        appendLinkToLinkOutputImpl(link, idMap);
        outputLock.unlock();
    }
    //! Write artificial sequences to disk
    void outputArtificialSequences(FastaCollection const & fastaCollection) {
        if (config_->outputArtificialSequences().string() != "") {
            std::unique_lock<std::mutex> outputLock(mutexOutput_);
            std::ofstream osArt(config_->outputArtificialSequences());
            if (osArt.good()) {
                for (auto&& rep : fastaCollection.collection()) {
                    rep.second.writeArtificialSequences(osArt);
                }
            } else {
                std::cerr << "[WARNING] -- Cannot write artificial sequences to file '" << config_->outputArtificialSequences() << "'" << std::endl;
            }
            outputLock.unlock();
        }
    }
    //! Write all Cubes to cube dict
    void outputCubes(Cubeset const & cubeset) {
        if (cubeDict_) {
            std::unique_lock<std::mutex> outputLock(mutexOutput_);
            auto const & idMap = *(cubeset.idMap());
            auto const & sequenceLenghts = *(cubeset.sequenceLengths());
            for (auto&& elem : cubeset.scoreToCube()) {
                auto score = elem.first;
                auto& cubes = elem.second;
                if (score >= config_->cubeScoreThreshold()) {
                    for (auto&& cubeptr : cubes) {
                        // create custom key to identify a cube
                        std::vector<JsonValue> keyVector;
                        for (auto&& td : cubeptr->tiledistance()) {
                            auto keyParts = std::array<std::string, 4>{idMap.queryGenomeName(td.genome()),
                                                                       idMap.querySequenceName(td.sequence()),
                                                                       std::to_string(td.distance()),
                                                                       std::to_string(td.reverse())};
                            keyVector.emplace_back(keyParts);
                        }
                        auto cubeKey = JsonValue(keyVector);
                        // collect links in cube and their counts
                        std::map<std::string, JsonValue> innerCubeMap;
                        auto& links = cubeset.cubeMap().at(cubeptr);
                        if (config_->cubeOutput() == 2) { // only on highest cube output level
                            for (auto&& link : links) {
                                std::vector<JsonValue> occurrences; // link id
                                for (auto&& occ : link.occurrence()) {
                                    occurrences.emplace_back(occ.toJsonValue(link.span(),
                                                                             idMap));
                                }
                                innerCubeMap.emplace(JsonValue(occurrences).value(),
                                                     cubeset.underlyingLinkset()->linkCount(link));
                            }
                        }
                        innerCubeMap.emplace("score", score);    // also include cube score
                        innerCubeMap.emplace("rawCount", links.size()); // and raw link count
                        // get region tuples
                        auto regionTuples = regionTupleExtraction(links,
                                                                  sequenceLenghts.at(cubeptr->tiledistance(0).sequence()),
                                                                  config_->tileSize(),
                                                                  1, 1, 1);
                        innerCubeMap.emplace("regionTuples", regionTuples);
                        cubeDict_->addValue(cubeKey.value(), innerCubeMap);
                    }
                }
            }
            outputLock.unlock();
        }
    }
    //! Write relevant cubes and their link counts to disk
    void outputRelevantCubes(PrefilterCubeset const & relevantCubeset) {
        std::unique_lock<std::mutex> outputLock(mutexOutput_);
        auto& idMap = *(relevantCubeset.idMap());
        for (auto&& elem : relevantCubeset.cubeMap()) {
            auto& cubeptr = elem.first;
            auto linkCount = elem.second;
            // create custom key to identify a cube
            std::vector<JsonValue> keyVector;
            for (auto&& td : cubeptr->tiledistance()) {
                auto keyParts = std::array<std::string, 4>{idMap.queryGenomeName(td.genome()),
                                                           idMap.querySequenceName(td.sequence()),
                                                           std::to_string(td.distance()),
                                                           std::to_string(td.reverse())};
                keyVector.emplace_back(keyParts);
            }
            auto cubeKey = JsonValue(keyVector);
            relevantCubeDict_->addValue(cubeKey.value(), linkCount);
        }
        outputLock.unlock();
    }
    template <typename LinkContainer>
    std::pair<size_t, size_t> outputLinks(LinkContainer const & links, IdentifierMapping const & idMap) {
        std::unique_lock<std::mutex> outputLock(mutexOutput_);
        size_t skipped = 0;
        size_t written = 0;
        for (auto&& link : links) {
            auto& occ0 = link.first();
            auto& occ1 = link.second();
            // only interested in matches between genome0 and genome1
            if (occ0.genome() > 1 || occ1.genome() > 1) {
                ++skipped;
                continue;
            }
            appendLinkToLinkOutputImpl(link, idMap);
            ++written;
        }
        outputLock.unlock();
        return std::pair<size_t, size_t>{written, skipped};
    }
    //! Append all Link s in a Linkset to link output
    template <typename LinksetType>
    std::pair<size_t, size_t> outputLinkset(LinksetType const & linkset, bool quiet = false) {
        std::unique_lock<std::mutex> outputLock(mutexOutput_);
        auto stats = outputLinksetImpl(linkset, quiet);
        outputLock.unlock();
        return stats;
    }
    //! Try to write all Links from linkset to output file, return true if succeeded
    template <typename LinksetType>
    bool tryOutputLinkset(LinksetType const & linkset, bool silent = true) {
        std::unique_lock<std::mutex> outputLock(mutexOutput_, std::defer_lock);
        if (outputLock.try_lock()) {
            outputLinksetImpl(linkset, silent);
            outputLock.unlock();
            return true;
        } else {
            return false;
        }
    }

private:
    //! Add sinlge link to output
    template <typename LinkType>
    void appendLinkToLinkOutputImpl(LinkType const & link,
                                    IdentifierMapping const & idMap) {
        auto span = link.span();
        std::vector<JsonValue> v;
        for (auto&& occ : link.occurrence()) { v.emplace_back(occ.toJsonValue(span, idMap)); }
        auto jsonVal = JsonValue(v);
        linkArray_ << jsonVal;
    }
    //! Append all Link s in a Linkset to link output
    template <typename LinksetType>
    std::pair<size_t, size_t> outputLinksetImpl(LinksetType const & linkset, bool quiet) {
        if (!quiet) { std::cout << "Writing matches to output file..." << std::endl; }
        size_t skipped = 0;
        size_t written = 0;
        for (auto&& elem : linkset.linkset()) {
            auto& link = elem.first;
            auto& occ0 = link.first();
            auto& occ1 = link.second();
            // only interested in matches between genome0 and genome1
            if (occ0.genome() > 1 || occ1.genome() > 1) {
                ++skipped;
                continue;
            }
            appendLinkToLinkOutputImpl(link, *(linkset.idMapping()));
            ++written;
        }
        if (!quiet) {
            std::cout << "[INFO] -- output -- Wrote " << written << " matches to output file" << std::endl;
            std::cout << "                    Skipped " << skipped << " matches that not include genome1 and 2" << std::endl << std::endl;
        }
        return std::pair<size_t, size_t>{written, skipped};
    }

    std::shared_ptr<Configuration const> config_;
    std::unique_ptr<JsonStreamDict> cubeDict_;
    std::unique_ptr<std::ofstream> cubeOs_;
    std::ofstream linkOs_;
    JsonStreamArray linkArray_;
    std::mutex mutexOutput_;
    std::unique_ptr<JsonStreamDict> relevantCubeDict_;
    std::unique_ptr<std::ofstream> relevantCubeOs_;
    std::unique_ptr<JsonStreamDict> runInfo_;
    std::unique_ptr<std::ofstream> runOs_;
};

#endif // OUTPUT_H
