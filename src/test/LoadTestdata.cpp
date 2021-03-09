#include "LoadTestdata.h"

LoadTestdataMetagraph::LoadTestdataMetagraph()
    : configFasta_{nullptr},
      configGraph_{nullptr},
      expectedGenomeNamesFasta_{},
      expectedGenomeNamesGraph_{},
      expectedSequenceLengthsExact_{std::make_shared<tsl::hopscotch_map<size_t, size_t>>()},
      expectedSequenceLengthsGraph_{std::make_shared<tsl::hopscotch_map<size_t, size_t>>()},
      expectedSequenceNamesFasta_{},
      expectedSequenceNamesGraph_{},
      fastas_{},
      graphAnnotationFile_{filepath_ / "testdataGraph.row.annodbg"},
      graphBinsize_{100},
      graphFile_{filepath_ / "testdataGraph.dbg"},
      idMapFasta_{nullptr},
      idMapGraph_{nullptr} {
    // fill fasta paths and expected genome and sequence names and lengths
    std::vector<size_t> seqlens;
    for (auto&& it : fs::directory_iterator(filepath_)) {
        auto entry = it.path();
        if (entry.has_extension() && (entry.extension() == ".fa" || entry.extension() == ".fasta")) {
            fastas_.emplace_back(entry);
            auto entryNoExt = entry;
            expectedGenomeNamesFasta_.emplace_back(entryNoExt.replace_extension().filename());
            expectedGenomeNamesGraph_.emplace_back(entry.filename());
            auto is = std::ifstream(entry);
            if (!is.is_open()) { throw std::runtime_error("[ERROR] -- LoadTestdataMetagraph -- could not open " + entry.string()); }
            std::string line{};
            while (is.good()) {
                std::getline(is, line);
                if (line.size() && line.at(0) == '>') {
                    expectedSequenceNamesFasta_.emplace_back(line.substr(1),
                                                             entryNoExt.filename());
                    expectedSequenceNamesGraph_.emplace_back(line.substr(1),
                                                             entry.filename());
                } else if (line.size()) {
                    seqlens.emplace_back(line.size());
                }
            }
        }
    }
    // make sure genome0.fa is first, i.e. getting genomeID 0
    std::sort(expectedGenomeNamesFasta_.begin(), expectedGenomeNamesFasta_.end());
    std::sort(expectedGenomeNamesGraph_.begin(), expectedGenomeNamesGraph_.end());
    if (seqlens.size() != expectedSequenceNamesFasta_.size()) { throw std::runtime_error("[ERROR] -- LoadTestdataMetagraph -- sequence lengths or expected sequences not correct"); }
    // create IDMaps and expectedSequenceLengths
    idMapFasta_ = std::make_shared<IdentifierMapping>(expectedGenomeNamesFasta_.at(0));
    idMapGraph_ = std::make_shared<IdentifierMapping>(expectedGenomeNamesGraph_.at(0));
    for (auto& genome : expectedGenomeNamesFasta_) {
        idMapFasta_->queryGenomeID(genome);
    }
    for (auto& genome : expectedGenomeNamesGraph_) {
        idMapGraph_->queryGenomeID(genome);
    }
    for (size_t i = 0; i < expectedSequenceNamesFasta_.size(); ++i) {
        auto& sequenceFasta = expectedSequenceNamesFasta_.at(i);
        auto& sequenceGraph = expectedSequenceNamesGraph_.at(i);
        auto len = seqlens.at(i);
        auto sid = idMapFasta_->querySequenceID(sequenceFasta.first, sequenceFasta.second);
        idMapGraph_->querySequenceID(sequenceGraph.first, sequenceGraph.second);
        (*expectedSequenceLengthsExact_)[sid] = len;
        (*expectedSequenceLengthsGraph_)[sid] = roundSeqLenToGraph(len) + 1;
    }
    // create configs for graph and fasta mode
    std::cout << "[INFO] -- LoadTestdataMetagraph -- Loading graph '" << graphFile_ << "' and annotation '" << graphAnnotationFile_ << "'" << std::endl;
    std::vector<std::string> fastastrings;
    for (auto& fa : fastas_) { fastastrings.emplace_back(fa.string()); }
    ConfigBuilder configBuilder{};
    configBuilder.set(Allvsall{true},
                      CubeLengthCutoff{0},
                      CubeScoreNormalizationParameter{3},
                      CubeScoreParameter{25},
                      CubeScoreParameterChunks{200},
                      CubeScoreThreshold{0},
                      GraphAnnotationFile{graphAnnotationFile_},
                      GraphFile{graphFile_},
                      Genome1{idMapGraph_->queryGenomeName(0)},
                      Genome2{idMapGraph_->queryGenomeName(1)},
                      Hasse(true),
                      InputFiles{fastastrings},
                      MaskCollectionPtr{std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight{18},
                                                                                         SpacedSeedMaskCollection::SeedSetSize{4})},
                      MetagraphInterfacePtr{std::make_shared<MetagraphInterface const>(graphFile_.string(),
                                                                                       graphAnnotationFile_.string())},
                      NThreads{1},
                      MatchLimit{ULLONG_MAX},
                      OptimalSeed{true},
                      PerformGeometricHashing{true},
                      Redmask{true},
                      Thinning{5},
                      TileSize{500});
    configGraph_ = configBuilder.makeConfig();
    // tweak fasta config
    configBuilder.set(Genome1{idMapFasta_->queryGenomeName(0)},
                      Genome2{idMapFasta_->queryGenomeName(1)});
    configFasta_ = configBuilder.makeConfig();
}



bool LoadTestdataMetagraph::checkExpectedGenomeNames(std::vector<std::string> const & thisGenomeNames,
                                                     std::vector<std::string> const & genomeNames) const {
    if (thisGenomeNames.size() != genomeNames.size()) { return false; }
    for (auto& genome : thisGenomeNames) {
        if (std::find(genomeNames.begin(), genomeNames.end(), genome) == genomeNames.end()) {
            std::cout << "[DEBUG] -- checkExpectedGenomeNames() -- '" << genome << "' not found in parameter genomeNames" << std::endl;
            return false;
        }
    }
    return true;
}



// check if all genome and sequence names are there in both, ids may differ
bool LoadTestdataMetagraph::checkIDMap(IdentifierMapping const & thisIDMap, IdentifierMapping const & idMap) const {
    if (thisIDMap.numGenomes() != idMap.numGenomes()) {
        std::cout << "[DEBUG] -- checkIDMap() -- Number of genomes differ in this (" << thisIDMap.numGenomes() << ") vs. parameter (" << idMap.numGenomes() << ")" << std::endl;
        return false;
    }
    if (thisIDMap.numSequences() != idMap.numSequences()) {
        std::cout << "[DEBUG] -- checkIDMap() -- Number of sequences differ in this (" << thisIDMap.numSequences() << ") vs. parameter (" << idMap.numSequences() << ")" << std::endl;
        return false;
    }
    for (auto& elem : thisIDMap.genomeIDToName()) {
        if (!idMap.genomeKnown(elem.second)) {
            std::cout << "[DEBUG] -- checkIDMap() -- Genome '" << elem.second << "' not known in parameter (1)" << std::endl;
            return false;
        }
    }
    for (auto& elem : thisIDMap.sequenceIDToTuple()) {
        auto genome = thisIDMap.queryGenomeName(elem.second.gid);
        if (!idMap.genomeKnown(genome)) {
            std::cout << "[DEBUG] -- checkIDMap() -- Genome '" << genome << "' not known in parameter (2)" << std::endl;
            return false;
        }
        if (!idMap.sequenceKnown(elem.second.sequence, genome)) {
            std::cout << "[DEBUG] -- checkIDMap() -- Sequence '" << elem.second.sequence << ", " << genome << "' not known in parameter" << std::endl;
            return false;
        }
    }
    for (auto& elem : idMap.genomeIDToName()) {
        if (!thisIDMap.genomeKnown(elem.second)) {
            std::cout << "[DEBUG] -- checkIDMap() -- Genome '" << elem.second << "' not known in this (1)" << std::endl;
            return false;
        }
    }
    for (auto& elem : idMap.sequenceIDToTuple()) {
        auto genome = idMap.queryGenomeName(elem.second.gid);
        if (!thisIDMap.genomeKnown(genome)) {
            std::cout << "[DEBUG] -- checkIDMap() -- Genome '" << genome << "' not known in this (2)" << std::endl;
            return false;
        }
        if (!thisIDMap.sequenceKnown(elem.second.sequence, genome)) {
            std::cout << "[DEBUG] -- checkIDMap() -- Sequence '" << elem.second.sequence << ", " << genome << "' not known in this" << std::endl;
            return false;
        }
    }
    return true;
}
// check if all sequence lengths are equal, ids may differ
bool LoadTestdataMetagraph::checkSequenceLengths(tsl::hopscotch_map<size_t, size_t> const & thisSeqlens,
                                                 IdentifierMapping const & thisIDMap,
                                                 tsl::hopscotch_map<size_t, size_t> const & seqlens,
                                                 IdentifierMapping const & idMap) const {
    if (seqlens.size() != thisSeqlens.size()) {
        std::cout << "[DEBUG] -- checkSequenceLengths() -- Sizes differ in this (" << thisSeqlens.size() << ") vs. parameter (" << seqlens.size() << ")" << std::endl;
        return false;
    }
    if (thisIDMap.numSequences() != idMap.numSequences()) {
        std::cout << "[DEBUG] -- checkSequenceLengths() -- Number of sequences differ in this idMap (" << thisIDMap.numSequences() << ") vs. parameter idMap (" << idMap.numSequences() << ")" << std::endl;
        return false;
    }
    for (auto& elem : seqlens) {
        auto& sequence = idMap.querySequenceTuple(elem.first);
        auto& genome = idMap.queryGenomeName(sequence.gid);
        if (!thisIDMap.genomeKnown(genome)) {
            std::cout << "[DEBUG] -- checkSequenceLengths() -- Genome '" << genome << "' not known in this idMap" << std::endl;
            return false;
        }
        if (!thisIDMap.sequenceKnown(sequence.sequence, genome)) {
            std::cout << "[DEBUG] -- checkSequenceLengths() -- Sequence '" << sequence.sequence << ", " << genome << "' not known in this idMap" << std::endl;
            return false;
        }
        auto sid = thisIDMap.querySequenceIDConst(sequence.sequence, genome);
        if (thisSeqlens.at(sid) != elem.second) {
            std::cout << "[DEBUG] -- checkSequenceLengths() -- Lengths for sequence '" << sequence.sequence << ", " << genome << "' differ in this (" << thisSeqlens.at(sid) << ") vs. parameter (" << elem.second << ")" << std::endl;
            return false;
        }
    }
    return true;
}
// load expected cubes
tsl::hopscotch_map<std::shared_ptr<Cube const>,
                   tsl::hopscotch_set<LinkPtr, LinkPtrHash>,
                   CubePtrHash, CubePtrEqual> LoadTestdataMetagraph::expectedCubes(fs::path const & file, IdentifierMapping const & idMap) const {
    std::cout << "[INFO] -- Loading file " << file << std::endl;
    tsl::hopscotch_map<std::shared_ptr<Cube const>,
                       tsl::hopscotch_set<LinkPtr, LinkPtrHash>,
                       CubePtrHash, CubePtrEqual> expectedCubeMap{};
    // load json
    std::ifstream is(file);
    json j;
    is >> j;
    // fill cubeMap
    for (auto&& elem : j.items()) {
        // create cube
        std::vector<Tiledistance> cubev;
        long long minDist = 0;
        for (auto&& occ : elem.value()["cube"]) {
            cubev.emplace_back(idMap.queryGenomeIDConst(occ[0].get<std::string>()),
                               idMap.querySequenceIDConst(occ[1].get<std::string>(), occ[0].get<std::string>()),
                               occ[2].get<long long>(),
                               occ[3].get<bool>());
            minDist = std::min(minDist, cubev.rbegin()->distance());
        }
        auto cubelink = Link();
        for (auto&& td : cubev) {
            cubelink.insertOccurrence(td.genome(), td.sequence(), (td.distance() + (-1*minDist)), td.reverse(), "A");
        }
        auto cube = std::make_shared<Cube const>(cubelink, configGraph_->tileSize());
        if (cube->tiledistance(0).distance() != 0) {
            std::cerr << "[ERROR] -- LoadTestdataMetagraph::expectedCubes -- cube '" << *cube << "' from '" << elem.value()["cube"] << "'" << std::endl;
            throw std::runtime_error("[ERROR] -- LoadTestdataMetagraph::expectedCubes -- Cube creation corrupted");
        }
        if (!(cubev == cube->tiledistance())) {
            std::cerr << "[ERROR] -- LoadTestdataMetagraph::expectedCubes -- cube vs. cubev: '" << cube << "' vs. '" << cubev << "'" << std::endl;
            throw std::runtime_error("[ERROR] -- LoadTestdataMetagraph::expectedCubes -- Cube creation corrupted");
        }
        // create links
        tsl::hopscotch_set<LinkPtr, LinkPtrHash> links;
        for (auto&& jlink : elem.value()["links"]) {
            auto link = Link();
            for (auto&& occ : jlink[0]) {
                link.insertOccurrence(idMap.queryGenomeIDConst(occ[0].get<std::string>()),
                                      idMap.querySequenceIDConst(occ[1].get<std::string>(), occ[0].get<std::string>()),
                                      occ[2].get<size_t>(),
                                      occ[3].get<bool>(),
                                      "A"); // kmer does not matter
            }
            link.extendSpanToRight(jlink[1].get<size_t>() - 1); // link initialized with span = 1
            links.emplace(LinkPtr(link));
        }
        // group links as links in expected cubes are not grouped in python but in SeedFinder pipeline
        /*std::vector<LinkPtr> sortedLinks(links.begin(), links.end());
        links.clear();
        std::sort(sortedLinks.begin(), sortedLinks.end());
        groupOverlappingSeeds(sortedLinks);
        links.insert(sortedLinks.begin(), sortedLinks.end());*/
        expectedCubeMap.emplace(cube, links);
    }
    return expectedCubeMap;
}
// load expected scores
tsl::hopscotch_map<double,
                   std::vector<std::shared_ptr<Cube const>>> LoadTestdataMetagraph::expectedCubeScores(fs::path const & file, IdentifierMapping const & idMap) const {
    std::cout << "[INFO] -- Loading file " << file << std::endl;
    tsl::hopscotch_map<double,
                       std::vector<std::shared_ptr<Cube const>>> expectedScoreMap{};
    // load json
    std::ifstream is(file);
    json j;
    is >> j;
    // fill cubeMap
    for (auto&& elem : j.items()) {
        // create cube
        std::vector<Tiledistance> cubev;
        long long minDist = 0;
        for (auto&& occ : elem.value()["cube"]) {
            cubev.emplace_back(idMap.queryGenomeIDConst(occ[0].get<std::string>()),
                               idMap.querySequenceIDConst(occ[1].get<std::string>(), occ[0].get<std::string>()),
                               occ[2].get<long long>(),
                               occ[3].get<bool>());
            minDist = std::min(minDist, cubev.rbegin()->distance());
        }
        auto cubelink = Link();
        for (auto&& td : cubev) {
            cubelink.insertOccurrence(td.genome(), td.sequence(), (td.distance() + (-1*minDist)), td.reverse(), "A");
        }
        auto cube = std::make_shared<Cube const>(cubelink, configGraph_->tileSize());
        if (cube->tiledistance(0).distance() != 0) {
            std::cerr << "[ERROR] -- LoadTestdataMetagraph::expectedCubeScores -- cube '" << *cube << "' from '" << elem.value()["cube"] << "'" << std::endl;
            throw std::runtime_error("[ERROR] -- LoadTestdataMetagraph::expectedCubeScores -- Cube creation corrupted");
        }
        if (!(cubev == cube->tiledistance())) {
            std::cerr << "[ERROR] -- LoadTestdataMetagraph::expectedCubeScores -- cube vs. cubev: '" << cube << "' vs. '" << cubev << "'" << std::endl;
            throw std::runtime_error("[ERROR] -- LoadTestdataMetagraph::expectedCubeScores -- Cube creation corrupted");
        }
        auto score = elem.value()["score"].get<double>();
        expectedScoreMap[score].emplace_back(cube);
    }
    return expectedScoreMap;
}
// load expected links
tsl::hopscotch_map<LinkPtr, size_t,
                   LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan> LoadTestdataMetagraph::expectedLinks(fs::path const & file, IdentifierMapping const & idMap) const {
    std::cout << "[INFO] -- Loading file " << file << std::endl;
    tsl::hopscotch_map<LinkPtr, size_t,
                       LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan> expectedLinkMap{};
    // load json
    std::ifstream is(file);
    json j;
    is >> j;
    // fill linkMap
    for (auto&& elem : j.items()) {
        auto link = Link();
        for (auto&& occ : elem.value()["link"]) {
            link.insertOccurrence(idMap.queryGenomeIDConst(occ[0].get<std::string>()),
                                  idMap.querySequenceIDConst(occ[1].get<std::string>(), occ[0].get<std::string>()),
                                  occ[2].get<size_t>(),
                                  occ[3].get<bool>(),
                                  "A"); // kmer does not matter
        }
        link.extendSpanToRight(elem.value()["span"].get<size_t>() - 1); // link initialized with span = 1
        expectedLinkMap.insert({LinkPtr(link), elem.value()["count"].get<size_t>()});
    }
    return expectedLinkMap;
}
// load expected region tuples
tsl::hopscotch_map<std::shared_ptr<Cube const>,
                   std::vector<std::array<size_t, 2>>,
                   CubePtrHash, CubePtrEqual> LoadTestdataMetagraph::expectedRegionTuples(fs::path const & file, IdentifierMapping const & idMap) const {
    std::cout << "[INFO] -- Loading file " << file << std::endl;
    tsl::hopscotch_map<std::shared_ptr<Cube const>,
                       std::vector<std::array<size_t, 2>>,
                       CubePtrHash, CubePtrEqual> expectedRegionTuples{};
    // load json
    std::ifstream is(file);
    json j;
    is >> j;
    // fill linkMap
    for (auto&& elem : j.items()) {
        // create cube
        std::vector<Tiledistance> cubev;
        long long minTile = 0;
        for (auto&& occ : elem.value()["cube"]) {
            auto tileID = occ[2].get<long long>();
            cubev.emplace_back(idMap.queryGenomeIDConst(occ[0].get<std::string>()),
                               idMap.querySequenceIDConst(occ[1].get<std::string>(), occ[0].get<std::string>()),
                               tileID,
                               occ[3].get<bool>());
            minTile = (tileID < minTile) ? tileID : minTile;  // add this to dummlink
        }
        auto dummylink = Link();  // fake link to create cube from
        for (auto& td : cubev) {
            dummylink.insertOccurrence(td.genome(),
                                       td.sequence(),
                                       (td.distance() - minTile) * configGraph_->tileSize(), // minTile is negative or zero
                                       td.reverse(),
                                       "A");
        }
        auto cube = std::make_shared<Cube const>(dummylink, configGraph_->tileSize());
        if (!(cubev == cube->tiledistance())) {
            std::cerr << "[ERROR] -- LoadTestdataMetagraph::expectedRegionTuples -- cube vs. cubev: '" << cube << "' vs. '" << cubev << "'" << std::endl;
            throw std::runtime_error("[ERROR] -- LoadTestdataMetagraph::expectedRegionTuples -- Cube creation corrupted");
        }
        // create tuples
        std::vector<std::array<size_t, 2>> tuples;
        for (auto&& jtuple : elem.value()["tuples"]) {
            tuples.emplace_back(std::array<size_t, 2>{jtuple[0].get<size_t>(),
                                                      jtuple[1].get<size_t>()});
        }
        expectedRegionTuples.emplace(cube, tuples);
    }
    return expectedRegionTuples;
}
// load expected seeds
tsl::hopscotch_map<TwoBitKmer<TwoBitKmerDataShort>,
                   std::vector<std::vector<KmerOccurrence>>,
                   TwoBitKmerHash<TwoBitKmerDataShort>> LoadTestdataMetagraph::expectedSeeds(fs::path const & file, IdentifierMapping const & idMap) const {
    std::cout << "[INFO] -- Loading file " << file << std::endl;
    tsl::hopscotch_map<TwoBitKmer<TwoBitKmerDataShort>,
                       std::vector<std::vector<KmerOccurrence>>,
                       TwoBitKmerHash<TwoBitKmerDataShort>> expectedSeedMap{};
    // load json
    std::ifstream is(file);
    json j;
    is >> j;
    // fill sorted seedMap
    for (auto&& elem : j.items()) {
        auto seed = TwoBitKmer<TwoBitKmerDataShort>(elem.key());
        std::vector<std::vector<KmerOccurrence>> occurrences{};
        for (auto&& perMaskOccs : elem.value()) {
            std::vector<KmerOccurrence> occs{};
            for (auto&& occ : perMaskOccs) {
                occs.emplace_back(idMap.queryGenomeIDConst(occ[0].get<std::string>()),
                                  idMap.querySequenceIDConst(occ[1].get<std::string>(), occ[0].get<std::string>()),
                                  occ[2].get<size_t>(),
                                  occ[3].get<bool>(),
                                  occ[4].get<std::string>());
            }
            std::sort(occs.begin(), occs.end());
            occurrences.emplace_back(occs);
        }
        expectedSeedMap.emplace(seed, occurrences);
    }
    return expectedSeedMap;
}
