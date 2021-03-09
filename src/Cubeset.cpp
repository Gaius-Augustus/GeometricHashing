#include "Cubeset.h"

Cubeset::Cubeset(std::shared_ptr<Linkset<LinkPtr,
                                         LinkPtrHashIgnoreSpan,
                                         LinkPtrEqualIgnoreSpan> const> linkset,
                 std::shared_ptr<tsl::hopscotch_map<size_t, size_t> const> sequenceLengths,
                 bool parallel)
    : cubeMap_{}, /*linkCount_{},*/linkFraction_{0}, linkset_{linkset},
      parallel_{parallel && linkset->config()->nThreads() > 1},
      sequenceLengths_{sequenceLengths},
      scoreToCube_{}, subcubeMap_{}, tilemap_{} {
    // determine if parallel or single thread (depends on allvsall parameter)
    size_t nThreads = (parallel_) ? config()->nThreads() : 1;
    // determine if running quiet
    bool quiet = !parallel_; // assume that when parallel execution is not allowed, this is called from multiple threads, thus no ouptu
    bool quietPb = quiet || config()->verbose() < 2;

    // Run Cubeset creation
    std::unique_ptr<MM::MemoryMonitor> mm;
    if (!quiet) {
        std::cout << "[INFO] -- Number of links: " << linkset_->size() << std::endl;
        mm = std::make_unique<MM::MemoryMonitor>();
        std::cout << "Memory usage before cube creation" << std::endl << *mm << std::endl;
    }

    std::unique_ptr<Timestep> ts;
    if (!quiet) { ts = std::make_unique<Timestep>("Creating Cubes"); }
    ProgressBar pb(linkset_->size(), quietPb);
    size_t i = 0;
    // iterate over all links and create cubes as needed
    for (auto&& linkToCount : linkset_->linkset()) {
        ++i; pb.update(i);
        auto&& linkptr = linkToCount.first;
        auto cube = std::make_shared<Cube>(*linkptr, config()->tileSize());
        // add new cubes to tilemap
        if (cubeMap_.find(cube) == cubeMap_.end()) {
            // only reference tiles are queried in subcube finding
            tilemap_[cube->tiledistance(0)].emplace_back(cube);
        }
        cubeMap_[cube].emplace(linkptr);
    }
    pb.finish();
    if (!quiet) {
        ts->endAndPrint();
        std::cout << cubeMap_.size() << " Cubes created" << std::endl << std::endl;
        std::cout << "Memory usage after cube creation" << std::endl << *mm << std::endl;
    }

    std::mutex mutex{};
    // find all subcubes of each cube (if hasse)
    auto fillSubcubes = [&mutex, this](SubcubeMapType & subcubeMap,
                                       typename CubeMapType::const_iterator it,
                                       typename CubeMapType::const_iterator end) {
        std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
        SubcubeMapType subcubeMapLocal;
        for (; it != end; ++it) {
            auto supercubes = getSupercubes(it->first);
            if (supercubes.size()) {
                for (auto&& supercube : supercubes) {
                    subcubeMapLocal[supercube].emplace(it->first);
                }
            }
        }
        lock.lock();
        for (auto&& elem : subcubeMapLocal) {
            subcubeMap[elem.first].insert(elem.second.begin(), elem.second.end());
        }
        lock.unlock();
    };
    if (config()->hasse() && numGenomes() > 2) {
        if (!quiet) {ts = std::make_unique<Timestep>("Finding Subcubes ("+std::to_string(nThreads)+" Threads)"); }
        executeParallel(cubeMap_, nThreads, fillSubcubes,
                        std::ref(subcubeMap_));
        if (!quiet) { ts->endAndPrint(); }
    }

    // compute cube scores
    //  compute scoring factor (number of links)/(number of possible links) = (nlinks)/prod(seqlens) [loop: avoid overflow]
    auto nlinks = linkset_->size();
    auto& idMap = *(linkset_->idMapping());
    std::vector<size_t> sequenceLenghtSums(idMap.numGenomes(), 0); // sum sequence lengths for each input genome
    for (auto&& elem : *sequenceLengths_) { sequenceLenghtSums.at(idMap.querySequenceTuple(elem.first).gid) += elem.second; }
    linkFraction_ = static_cast<long double>(nlinks);
    for (auto lengthsum : sequenceLenghtSums) { linkFraction_ /= static_cast<long double>(lengthsum); }
    if (!quiet) {
        std::cout << std::endl << "[INFO] -- number of links: " << linkset_->size() << std::endl;
        std::cout <<              "[INFO] -- linkFraction_: " << linkFraction_ << std::endl << std::endl;
    }
    auto computeCubeScores = [&mutex, this](ScoreMapType & scoreToCube,
                                            typename CubeMapType::const_iterator it,
                                            typename CubeMapType::const_iterator end) {
        std::unique_lock<std::mutex> lock(mutex, std::defer_lock);
        ScoreMapType scoreToCubeLocal;
        for (; it != end; ++it) {
            if (config()->oldCubeScore()) {
                scoreToCubeLocal[computeCubeScoreOld(it->first)].emplace_back(it->first);
            } else {
                scoreToCubeLocal[computeCubeScore(it->first)].emplace_back(it->first);
            }
        }
        lock.lock();
        for (auto&& elem : scoreToCubeLocal) {
            scoreToCube[elem.first].insert(scoreToCube[elem.first].end(),
                                           scoreToCubeLocal[elem.first].begin(),
                                           scoreToCubeLocal[elem.first].end());
        }
        lock.unlock();
    };
    if (!quiet) {ts = std::make_unique<Timestep>("Computing Cube scores ("+std::to_string(nThreads)+" Threads)"); }
    executeParallel(cubeMap_, nThreads, computeCubeScores,
                    std::ref(scoreToCube_));
    if (!quiet) { ts->endAndPrint(); }

    if (!quiet) { std::cout << "Memory usage after Cubeset creation" << std::endl << *mm << std::endl; }
}



double Cubeset::computeCubeScore(std::shared_ptr<Cube const> cube) const {
    size_t tilesize = config()->tileSize();
    // if no tilesize, set tilesize to biggest sequence such that entire cube "fits" in single tile
    if (tilesize == 0) {
        for (auto&& td : cube->tiledistance()) {
            tilesize = std::max(tilesize, sequenceLengths_->at(td.sequence()));
        }
    }
    auto chunksize = config()->cubeScoreParameterChunks();
    auto p = config()->cubeScoreNormalizationParameter();

    // if tilesize mod nsubtilesPerSide > 0, subtilewidth * subtilesize > tilesize, adjust these values below
    auto dimExp = static_cast<long double>(cube->dimensionality()-1);
    auto nsubtilesPerSide = (dimExp > 1)
            ? static_cast<size_t>(std::floor(std::pow(nsubtiles(), 1.L/dimExp)))
            : nsubtiles();
    auto subtilewidth = tilesize / nsubtilesPerSide; // both size_t, thus result in floor(tilesize/nsubtilesPerSide)
    nsubtilesPerSide += (tilesize % nsubtilesPerSide == 0) ? 0 : 1; // one subtile more per side if subtilewidth*nsubtilesPerSide < tilesize
    auto nsubtiles = static_cast<size_t>(std::pow(nsubtilesPerSide, dimExp));

    auto cubeGeo = cube->geoProperties(tilesize, chunksize, *sequenceLengths_);
    auto nchunks =  (chunksize) ? cubeGeo.maxChunk - cubeGeo.minChunk + 1 : 1ull;
    auto chunklen = (chunksize) ? chunksize : std::max((size_t)1, cubeGeo.length); // length may be zero
    auto chunkvol = chunklen * static_cast<size_t>(std::pow(subtilewidth, dimExp));
    auto lambda = static_cast<long double>(chunkvol) * linkFraction_;

    tsl::hopscotch_map<Cube, std::vector<LinkPtr>, CubeHash> subcubes;
    for (auto&& link : cubeMap_.at(cube)) {
        Cube subcube(*link, subtilewidth);
        subcubes[subcube].emplace_back(link);
    }

    size_t pnorm = 0;
    for (auto&& elem : subcubes) {
        tsl::hopscotch_map<size_t, size_t> chunkCounts;
        for (auto&& link : elem.second) {
            auto chunkID = link.chunkID(chunksize);
            if (chunkCounts.find(chunkID) == chunkCounts.end()) { chunkCounts[chunkID] = 0; }
            chunkCounts[chunkID] += 1;
        }
        for (auto&& count : chunkCounts) { pnorm += static_cast<size_t>(std::pow(count.second, p)); }
    }

    auto normalization = lambda * std::pow(static_cast<long double>(nsubtiles*nchunks), 1.L/static_cast<long double>(p));
    auto score = std::pow(pnorm, 1.L/static_cast<long double>(p)) / normalization;
    if (score < 0.) {
        std::cerr << "[DEBUG] -- cube " << *cube << std::endl;
        std::cerr << "           tilesize " << tilesize << std::endl;
        std::cerr << "           chunksize " << chunksize << std::endl;
        std::cerr << "           p " << p << std::endl;
        std::cerr << "           nsubtilesPerSide " << nsubtilesPerSide << std::endl;
        std::cerr << "           subtilesize " << subtilewidth << std::endl;
        std::cerr << "           nsubtiles " << nsubtiles << std::endl;
        std::cerr << "           cubeGeo " << cubeGeo.length << ", " << cubeGeo.minChunk << ", " << cubeGeo.maxChunk << std::endl;
        std::cerr << "           nchunks " << nchunks << std::endl;
        std::cerr << "           chunklen " << chunklen << std::endl;
        std::cerr << "           chunkvol " << chunkvol << std::endl;
        std::cerr << "           seqlens " << *sequenceLengths_ << std::endl;
        std::cerr << "           linkFraction " << linkFraction_ << std::endl;
        std::cerr << "           lambda " << lambda << std::endl;
        std::cerr << "           pnorm " << pnorm << std::endl;
        std::cerr << "           normalization " << normalization << std::endl;
        std::cerr << "           score " << score << std::endl;
    }
    return static_cast<double>(score);
}



double Cubeset::computeCubeScoreOld(std::shared_ptr<Cube const> cube) const {
    (void)cube;
    throw std::runtime_error("[ERROR] -- Cubeset::computeCubeScoreOld -- deprecated");
    return -1.;
}



tsl::hopscotch_set<LinkPtr,
                   LinkPtrHash> Cubeset::linksIncludingSubcubes(std::shared_ptr<Cube const> const & cube) const {
    tsl::hopscotch_set<LinkPtr, LinkPtrHash> links{cubeMap_.at(cube).begin(), cubeMap_.at(cube).end()};
    if (subcubeMap_.find(cube) != subcubeMap_.end()) {
        for (auto&& subcube : subcubeMap_.at(cube)) {
            links.insert(cubeMap_.at(subcube).begin(), cubeMap_.at(subcube).end());
        }
    }
    return links;
}
