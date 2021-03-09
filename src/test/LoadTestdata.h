#ifndef LOADTESTDATA_H
#define LOADTESTDATA_H

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"
#include "../Cubeset.h"
#include "../ExtractSeeds.h"
#include "../FastaCollection.h"
#include "../FastaRepresentation.h"
#include "../IdentifierMapping.h"
#include "../Linkset.h"
#include "../regionTupleExtraction.h"
#include "../SeedMap.h"
#include "../TwoBitKmer.h"
#include "ConfigurationGenerator.h"
#include "GeometricHashingTestdata.h"

using json = nlohmann::json;



inline auto testfilesBasepath() {
    return fs::path(SEEDFINDING_TESTFILEPATH);
}


class LoadTestdata {
public:
    static struct GeometricHashingData{} geometricHashingData; // tag dispatch

    LoadTestdata()
        : cbuilder_{},
          fastaCollection_{std::make_shared<FastaCollection>()},
          genome0_{"hg38_orthologs"},
          genome1_{"mm10_orthologs"},
          idMap_{std::make_shared<IdentifierMapping>(genome0_)},    // id 0
          inputFiles_{filepath_ + "/hetGla2_orthologs.fa",
                      filepath_ + "/hg38_orthologs.fa",
                      filepath_ + "/macFas5_orthologs.fa",
                      filepath_ + "/mm10_orthologs.fa"},
          sequenceLengths_{std::make_shared<tsl::hopscotch_map<size_t, size_t>>()} {

        idMap_->queryGenomeID(genome1_);    // id 1
        for (auto&& file : inputFiles_) {
            fastaCollection_->emplace(FastaRepresentation::genomeFromFilename(file), file);
        }
        fastaCollection_->populateIdentifierMappingFromFastaCollection(*idMap_);
        for (auto&& sid : idMap_->sequenceIDToTuple()) {
            sequenceLengths_->emplace(sid.first, fastaCollection_->sequenceLength(sid.first, *idMap_));
        }
        cbuilder_.set(Genome1{genome0_},
                      Genome2{genome1_},
                      InputFiles{inputFiles_});
        std::cout << "Current working dir: " << std::filesystem::current_path() << std::endl;
    }
    LoadTestdata(GeometricHashingData)
        : cbuilder_{},
          fastaCollection_{std::make_shared<FastaCollection>()},
          genome0_{kTestSpecies[0]},
          genome1_{kTestSpecies[1]},
          idMap_{std::make_shared<IdentifierMapping>(genome0_)},    // id 0
          inputFiles_{filepath_ + "/geometricHashing/species1.fasta",
                      filepath_ + "/geometricHashing/species2.fasta",
                      filepath_ + "/geometricHashing/species3.fasta",
                      filepath_ + "/geometricHashing/species4.fasta",
                      filepath_ + "/geometricHashing/species5.fasta"},
          sequenceLengths_{std::make_shared<tsl::hopscotch_map<size_t, size_t>>()}
    {
        idMap_->queryGenomeID(genome1_);    // id 1
        idMap_->queryGenomeID(kTestSpecies[2]);    // id 2
        idMap_->queryGenomeID(kTestSpecies[3]);    // id 3
        idMap_->queryGenomeID(kTestSpecies[4]);    // id 4
        idMap_->querySequenceID(kTestSequences[0], kTestSpecies[0]);    // id 0
        idMap_->querySequenceID(kTestSequences[1], kTestSpecies[1]);    // id 1
        idMap_->querySequenceID(kTestSequences[2], kTestSpecies[1]);    // id 2
        idMap_->querySequenceID(kTestSequences[3], kTestSpecies[2]);    // id 3
        idMap_->querySequenceID(kTestSequences[4], kTestSpecies[2]);    // id 4
        idMap_->querySequenceID(kTestSequences[5], kTestSpecies[3]);    // id 5
        idMap_->querySequenceID(kTestSequences[6], kTestSpecies[3]);    // id 6
        idMap_->querySequenceID(kTestSequences[7], kTestSpecies[3]);    // id 7
        idMap_->querySequenceID(kTestSequences[8], kTestSpecies[4]);    // id 8
        idMap_->querySequenceID(kTestSequences[9], kTestSpecies[4]);    // id 9

        for (auto&& file : inputFiles_) {
            fastaCollection_->emplace(FastaRepresentation::genomeFromFilename(file), file);
        }
        for (auto&& sid : idMap_->sequenceIDToTuple()) {
            sequenceLengths_->emplace(sid.first, fastaCollection_->sequenceLength(sid.first, *idMap_));
        }
        cbuilder_.set(Genome1{genome0_},
                      Genome2{genome1_},
                      InputFiles{inputFiles_},
                      PerformGeometricHashing{true});
        std::cout << "Current working dir: " << std::filesystem::current_path() << std::endl;
    }

    auto & cbuilder() { return cbuilder_; }
    auto fastaCollection() const { return fastaCollection_; }
    auto const & filepath() const { return filepath_; }
    auto const & genome0() const { return genome0_; }
    auto const & genome1() const { return genome1_; }
    auto idMap() const { return idMap_; }
    auto const & inputFiles() const { return inputFiles_; }
    auto sequenceLengths() const { return sequenceLengths_; }

private:
    ConfigBuilder cbuilder_;
    std::shared_ptr<FastaCollection> fastaCollection_;
    std::string const filepath_ = SEEDFINDING_TESTFILEPATH;
    std::string genome0_;
    std::string genome1_;
    std::shared_ptr<IdentifierMapping> idMap_;
    std::vector<std::string> inputFiles_;
    std::shared_ptr<tsl::hopscotch_map<size_t, size_t>> sequenceLengths_;
};



inline bool discardKmer(std::string const & kmer, size_t thinning) {
    auto hash = std::hash<std::string>{}(kmer);
    return ((hash % thinning) == 1);
}



class LoadTestdataMetagraph {
public:
    LoadTestdataMetagraph();

    bool checkExpectedGenomeNamesFasta(std::vector<std::string> const & genomeNames) {
        return checkExpectedGenomeNames(expectedGenomeNamesFasta_, genomeNames);
    }
    bool checkExpectedGenomeNamesGraph(std::vector<std::string> const & genomeNames) {
        return checkExpectedGenomeNames(expectedGenomeNamesGraph_, genomeNames);
    }
    bool checkIDMapFasta(IdentifierMapping const & idMap) const { return checkIDMap(*idMapFasta_, idMap); }
    bool checkIDMapGraph(IdentifierMapping const & idMap) const { return checkIDMap(*idMapGraph_, idMap); }
    bool checkSequenceLengthsFasta(tsl::hopscotch_map<size_t, size_t> const & seqlens, IdentifierMapping const & idMap) const {
        return checkSequenceLengths(*expectedSequenceLengthsExact_, *idMapFasta_, seqlens, idMap);
    }
    bool checkSequenceLengthsGraph(tsl::hopscotch_map<size_t, size_t> const & seqlens, IdentifierMapping const & idMap) const {
        return checkSequenceLengths(*expectedSequenceLengthsGraph_, *idMapGraph_, seqlens, idMap);
    }
    auto configFasta() const { return configFasta_; }
    auto configGraph() const { return configGraph_; }
    auto const & expectedGenomeNamesFasta() const { return expectedGenomeNamesFasta_; }
    auto const & expectedGenomeNamesGraph() const { return expectedGenomeNamesGraph_; }
    // needs genome names with ".fa"
    auto expectedCubesFasta(IdentifierMapping const & graphIDMap) { return expectedCubes(filepath_ / "expectedCubes.json", graphIDMap); }
    auto expectedCubesGraph(IdentifierMapping const & graphIDMap) { return expectedCubes(filepath_ / "expectedCubesGraph.json", graphIDMap); }
    auto expectedCubesGroupedFasta(IdentifierMapping const & graphIDMap) { return expectedCubes(filepath_ / "expectedCubesGrouped.json", graphIDMap); }
    auto expectedCubesGroupedGraph(IdentifierMapping const & graphIDMap) { return expectedCubes(filepath_ / "expectedCubesGroupedGraph.json", graphIDMap); }
    auto expectedCubeScoresFasta(IdentifierMapping const & graphIDMap) { return expectedCubeScores(filepath_ / "expectedCubeScores.json", graphIDMap); }
    auto expectedCubeScoresGraph(IdentifierMapping const & graphIDMap) { return expectedCubeScores(filepath_ / "expectedCubeScoresGraph.json", graphIDMap); }
    auto expectedCubesWithSubcubeLinksFasta(IdentifierMapping const & graphIDMap) { return expectedCubes(filepath_ / "expectedCubesWithSubcubeLinks.json", graphIDMap); }
    auto expectedCubesWithSubcubeLinksGraph(IdentifierMapping const & graphIDMap) { return expectedCubes(filepath_ / "expectedCubesWithSubcubeLinksGraph.json", graphIDMap); }
    auto expectedCubesGroupedWithSubcubeLinksFasta(IdentifierMapping const & graphIDMap) { return expectedCubes(filepath_ / "expectedCubesGroupedWithSubcubeLinks.json", graphIDMap); }
    auto expectedCubesGroupedWithSubcubeLinksGraph(IdentifierMapping const & graphIDMap) { return expectedCubes(filepath_ / "expectedCubesGroupedWithSubcubeLinksGraph.json", graphIDMap); }
    // needs genome names with ".fa"
    auto expectedLinksFasta(IdentifierMapping const & graphIDMap) { return expectedLinks(filepath_ / "expectedLinks.json", graphIDMap); }
    auto expectedLinksGraph(IdentifierMapping const & graphIDMap) { return expectedLinks(filepath_ / "expectedLinksGraph.json", graphIDMap); }
    // needs genome names with ".fa"
    auto expectedRegionTuplesFasta(IdentifierMapping const & graphIDMap) { return expectedRegionTuples(filepath_ / "expectedRegionTuples.json", graphIDMap); }
    auto expectedRegionTuplesGraph(IdentifierMapping const & graphIDMap) { return expectedRegionTuples(filepath_ / "expectedRegionTuplesGraph.json", graphIDMap); }
    // needs genome names with ".fa"
    auto expectedSeedsFasta(IdentifierMapping const & graphIDMap) { return expectedSeeds(filepath_ / "expectedSeeds.json", graphIDMap); }
    auto expectedSeedsGraph(IdentifierMapping const & graphIDMap) { return expectedSeeds(filepath_ / "expectedSeedsGraph.json", graphIDMap); }
    std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> expectedSequenceLengthsExact() const { return expectedSequenceLengthsExact_; }
    std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> expectedSequenceLengthsGraph() const { return expectedSequenceLengthsGraph_; }
    auto const & expectedSequenceNamesFasta() const { return expectedSequenceNamesFasta_; }
    auto const & expectedSequenceNamesGraph() const { return expectedSequenceNamesGraph_; }
    auto const & fastas() const { return fastas_; }
    auto filepath() const { return filepath_; }
    auto graphAnnotationFile() const { return graphAnnotationFile_; }
    std::shared_ptr<IdentifierMapping const> idMapFasta() const { return idMapFasta_; }
    std::shared_ptr<IdentifierMapping const> idMapGraph() const { return idMapGraph_; }
    auto masks() const { return configGraph_->maskCollection(); }
    std::shared_ptr<MetagraphInterface const> metagraphInterface() const {
        return configGraph_->metagraphInterface();
    }
    size_t roundPositionToGraph(size_t pos) const {
        auto graphBinsize = static_cast<double>(graphBinsize_);
        return static_cast<size_t>(std::floor(static_cast<double>(pos)/graphBinsize)*graphBinsize);
    }
    size_t roundSeqLenToGraph(size_t trueLen) {
        if (trueLen == 0) {
            std::cerr << "[WARNING] -- LoadTestdataMetagraph::roundSeqLen -- true length is 0, returning 0" << std::endl;
            return 0;
        }
        return roundPositionToGraph(trueLen-1);
    }
    auto runSeedExtractionFasta(std::shared_ptr<IdentifierMapping> idMap,
                                std::shared_ptr<tsl::hopscotch_map<size_t, size_t>> sequenceLengths) const {
        auto fastaCollection = std::make_shared<FastaCollection>(configFasta_);
        fastaCollection->populateIdentifierMappingFromFastaCollection(*idMap);
        fastaCollection->fillSequenceLengths(*sequenceLengths, *idMap);
        auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(configFasta_, idMap);
        seedMap->extractSeeds(fastaCollection);
        return seedMap;
    }
    auto runSeedExtractionFasta() { return runSeedExtractionFasta(idMapFasta_, expectedSequenceLengthsExact_); }
    auto runSeedExtractionGraph(std::shared_ptr<IdentifierMapping> idMap,
                                std::shared_ptr<tsl::hopscotch_map<size_t, size_t>> sequenceLengths) const {
        auto seedMap = std::make_shared<SeedMap<TwoBitKmerDataShort>>(configGraph_, idMap);
        auto tsSingle = mabl3::Timestep("Single threaded seed extraction");
        //ExtractFromGraph<TwoBitKmerDataMedium, TwoBitKmerDataShort>(idMap,
        //                                                            sequenceLengths,
        //                                                            seedMap);
        auto callback = [this,
                         &seedMap,
                         &sequenceLengths,
                         &idMap](std::string const & kmerstr,
                                 std::vector<MetagraphInterface::NodeAnnotation> const & occurrences){
            for (auto&& annot : occurrences) {
                auto sid = idMap->querySequenceID(annot.sequence, annot.genome);
                if (sequenceLengths->find(sid) == sequenceLengths->end()) {
                    (*sequenceLengths)[sid] = annot.bin_idx;
                } else if (sequenceLengths->at(sid) < annot.bin_idx) {
                    (*sequenceLengths)[sid] = annot.bin_idx;
                }
            }
            // for now, deliberately (potentially) miss the last few bases of each sequence
            //   if k > max_span (see kmerstr.substr() below:
            //                      If this is not done, seedFinding will crash if max_span
            //                      induces a smaller TwoBitKmerData type than k!)
            // For fasta inputs and masks with span < max_span, last bases of each sequence
            //   are already lost and will be lost here as well
            // Further, shifting the masks over the k-mer could not guarantee the correct
            //   bin (tile) as this is determined only by the first base of the k-mer
            if (discardKmer(kmerstr, configGraph_->thinning())) { return; }
            //auto validKmer = std::make_shared<bool>();
            //auto kmer = TwoBitKmer<TwoBitKmerDataMedium>(kmerstr.substr(0, seedMap->config()->span()),
            //                                             validKmer);

            auto wrapper = [&seedMap](TwoBitKmer<TwoBitKmerDataShort> const & seed, KmerOccurrence const & occurrence, size_t maskIndex) {
                seedMap->addSeed(seed, occurrence, maskIndex);
            };

            ExtractSeeds<TwoBitKmerDataShort> extractor{configGraph_->maskCollection(), idMap, configGraph_};
            //if (*validKmer) {
                for (auto&& occ : occurrences) {
                    auto gid = idMap->queryGenomeIDConst(occ.genome);
                    auto sid = idMap->querySequenceIDConst(occ.sequence, occ.genome);
                    KmerOccurrence occurrence{(uint8_t) gid,
                                              (uint32_t) sid,
                                              occ.bin_idx,
                                              occ.reverse_strand,
                                              kmerstr};
                    extractor.createSeed<KmerOccurrence>(kmerstr.substr(0, seedMap->config()->span()), occurrence, wrapper);
                }
            //}
        };
        seedMap->config()->metagraphInterface()->iterateNodes(callback,0,0,false,false);
        tsSingle.endAndPrint();
        return seedMap;
    }
    auto runSeedExtractionGraph() { return runSeedExtractionGraph(idMapGraph_, expectedSequenceLengthsGraph_); }
    std::shared_ptr<Linkset<LinkPtr,
                            LinkPtrHashIgnoreSpan,
                            LinkPtrEqualIgnoreSpan> const> runLinksetCreation(SeedMap<TwoBitKmerDataShort> const & seedMap,
                                                                              std::shared_ptr<IdentifierMapping const> idMap) const {
        auto linkset = std::make_shared<Linkset<LinkPtr,
                                                LinkPtrHashIgnoreSpan,
                                                LinkPtrEqualIgnoreSpan>>(configGraph_,
                                                                         idMap); // configs don't differ here
        linkset->createLinks(seedMap);
        std::cout << "[INFO] -- LoadTestdataMetagraph::runLinksetCreation() -- Created " << linkset->size() << " Links" << std::endl;
        linkset->printStatistics(std::cout);
        return linkset;
    }
    auto runLinksetCreationFasta() {
        return runLinksetCreation(*runSeedExtractionFasta(), idMapFasta_);
    }
    auto runLinksetCreationGraph() {
        return runLinksetCreation(*runSeedExtractionGraph(), idMapGraph_);
    }
    auto runCubesetCreation(std::shared_ptr<Linkset<LinkPtr,
                                                    LinkPtrHashIgnoreSpan,
                                                    LinkPtrEqualIgnoreSpan> const> linkset,
                            std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> sequenceLengths) const {
        Cubeset cubeset(linkset, sequenceLengths);
        std::cout << "[INFO] -- LoadTestdataMetagraph::runCubesetCreation() -- Created " << cubeset.cubeMap().size() << " Cubes" << std::endl;
        return cubeset;
    }
    auto runCubesetCreationFasta() {
        return runCubesetCreation(runLinksetCreationFasta(), expectedSequenceLengthsExact_);
    }
    auto runCubesetCreationGraph() {
        return runCubesetCreation(runLinksetCreationGraph(), expectedSequenceLengthsGraph_);
    }
    auto runRegionTupleExtraction(Cubeset cubeset,
                                  tsl::hopscotch_map<size_t, size_t> const & sequenceLengths) const {
        cubeset.groupLinksInCubes();
        tsl::hopscotch_map<std::shared_ptr<Cube const>,
                           std::vector<std::array<size_t, 2>>,
                           CubePtrHash, CubePtrEqual> tuplesMap;
        for (auto& cube : cubeset.cubeMap()) {
            tuplesMap.emplace(cube.first, regionTupleExtraction(cube.second,
                                                                sequenceLengths.at(cube.first->tiledistance(0).sequence()),
                                                                configGraph_->tileSize(),
                                                                1, 1, 1));
        }
        return tuplesMap;
    }
    auto runRegionTupleExtractionFasta() {
        return runRegionTupleExtraction(runCubesetCreationFasta(), *expectedSequenceLengthsExact_);
    }
    auto runRegionTupleExtractionGraph() {
        return runRegionTupleExtraction(runCubesetCreationGraph(), *expectedSequenceLengthsGraph_);
    }
    static auto sortedSeedMap(SeedMap<TwoBitKmerDataShort> const & seedMap) {
        tsl::hopscotch_map<TwoBitKmer<TwoBitKmerDataShort>,
                           std::vector<std::vector<KmerOccurrence>>,
                           TwoBitKmerHash<TwoBitKmerDataShort>> sortedSeedMap{};
        for (auto& seed : seedMap.seedMap()) {
            std::vector<std::vector<KmerOccurrence>> sortedOuter;
            for (auto& inner : seed.second) {
                auto vector = inner;
                std::sort(vector.begin(), vector.end());
                sortedOuter.emplace_back(vector);
            }
            sortedSeedMap.emplace(seed.first, sortedOuter);
        }
        return sortedSeedMap;
    }
private:
    bool checkExpectedGenomeNames(std::vector<std::string> const & thisGenomeNames,
                                  std::vector<std::string> const & genomeNames) const;
    bool checkIDMap(IdentifierMapping const & thisIDMap, IdentifierMapping const & idMap) const;
    bool checkSequenceLengths(tsl::hopscotch_map<size_t, size_t> const & thisSeqlens,
                              IdentifierMapping const & thisIDMap,
                              tsl::hopscotch_map<size_t, size_t> const & seqlens,
                              IdentifierMapping const & idMap) const;
    tsl::hopscotch_map<std::shared_ptr<Cube const>,
                       tsl::hopscotch_set<LinkPtr, LinkPtrHash>,
                       CubePtrHash, CubePtrEqual> expectedCubes(fs::path const & file, IdentifierMapping const & idMap) const;
    tsl::hopscotch_map<double,
                       std::vector<std::shared_ptr<Cube const>>> expectedCubeScores(fs::path const & file, IdentifierMapping const & idMap) const;
    tsl::hopscotch_map<LinkPtr,
                       size_t,
                       LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan> expectedLinks(fs::path const & file, IdentifierMapping const & idMap) const;
    tsl::hopscotch_map<std::shared_ptr<Cube const>,
                       std::vector<std::array<size_t, 2>>,
                       CubePtrHash, CubePtrEqual> expectedRegionTuples(fs::path const & file, IdentifierMapping const & idMap) const;
    tsl::hopscotch_map<TwoBitKmer<TwoBitKmerDataShort>,
                       std::vector<std::vector<KmerOccurrence>>,
                       TwoBitKmerHash<TwoBitKmerDataShort>> expectedSeeds(fs::path const & file, IdentifierMapping const & idMap) const;

    std::shared_ptr<Configuration const> configFasta_;
    std::shared_ptr<Configuration const> configGraph_;
    std::vector<std::string> expectedGenomeNamesFasta_;
    std::vector<std::string> expectedGenomeNamesGraph_;
    std::shared_ptr<tsl::hopscotch_map<size_t, size_t>> expectedSequenceLengthsExact_;
    std::shared_ptr<tsl::hopscotch_map<size_t, size_t>> expectedSequenceLengthsGraph_;
    std::vector<std::pair<std::string, std::string>> expectedSequenceNamesFasta_;
    std::vector<std::pair<std::string, std::string>> expectedSequenceNamesGraph_;
    std::vector<fs::path> fastas_;
    fs::path const filepath_ = /*fs::path(SEEDFINDING_TESTFILEPATH)*/ testfilesBasepath() / "metagraph";
    fs::path graphAnnotationFile_;
    size_t graphBinsize_;
    fs::path graphFile_;
    std::shared_ptr<IdentifierMapping> idMapFasta_;
    std::shared_ptr<IdentifierMapping> idMapGraph_;
};



template <typename T>
inline bool hopscotchEqual(T const & h1, T const & h2) {
    auto ret = true;
    for (auto& elem : h1) {
        if (h2.find(elem.first) == h2.end()) {
            std::cout << "'" << elem.first << "' not found in h2" << std::endl;
            ret = false;//return false;
        }
        if (h2.at(elem.first) != elem.second) {
            std::cout << "For '" << elem.first << "', " << h2.at(elem.first) << " (h2) != " << elem.second << " (h1)" << std::endl;
            ret = false;//return false;
        }
    }
    for (auto& elem : h2) {
        if (h1.find(elem.first) == h1.end()) {
            std::cout << "'" << elem.first << "' not found in h1" << std::endl;
            ret = false;//return false;
        }
        if (h1.at(elem.first) != elem.second) {
            std::cout << "For '" << elem.first << "', " << h1.at(elem.first) << " (h1) != " << elem.second << " (h2)" << std::endl;
            ret = false;//return false;
        }
    }
    return ret;//return true;
}

template <typename T>
inline bool hopscotchSetEqual(T const & h1, T const & h2) {
    auto ret = true;
    for (auto& elem : h1) {
        if (h2.find(elem) == h2.end()) {
            std::cout << "'" << elem << "' not found in h2" << std::endl;
            ret = false;
        }
    }
    for (auto& elem : h2) {
        if (h1.find(elem) == h1.end()) {
            std::cout << "'" << elem << "' not found in h1" << std::endl;
            ret = false;
        }
    }
    return ret;
}

#endif // LOADTESTDATA_H
