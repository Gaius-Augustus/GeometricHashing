#include "Cubeset.h"

Cubeset::Cubeset(std::shared_ptr<Linkset const> linkset,
                 std::shared_ptr<FastaCollection const> fastaCollection,
                 std::shared_ptr<Configuration const> config)
    : config_{config}, cubeAreaCutoff_{config->cubeAreaCutoff()}, cubeMap_{},
      cubeScoreParameter_{config->cubeScoreParameter()},
      eta_{config->cubeScoreNormalizationParameter()},
      fastaCollection_{fastaCollection}, linkset_{linkset},
      numGenomes_{linkset->numGenomes()}, scoreToCube_{} {
    auto mm = MM::MemoryMonitor();
    std::cout << "Memory usage before cube creation" << std::endl << mm << std::endl;

    Timestep tsCube("Creating Cubes");
    ProgressBar pb(linkset_->numLinks(), config_->quiet());
    size_t i = 0;
    // iterate over all links and create cubes as needed
    for (auto&& linkToCount : linkset_->linkset()) {
        ++i; pb.update(i);
        auto&& linkptr = linkToCount.first;
        auto cube = std::make_shared<Cube>(linkptr, config->tileSize());
        cubeMap_[cube].emplace(linkptr);
    }
    std::cout << std::endl;
    tsCube.endAndPrint(Timestep::minutes);
    std::cout << cubeMap_.size() << " Cubes created" << std::endl << std::endl;
    std::cout << "Memory usage after cube creation" << std::endl << mm << std::endl;

    Timestep tsScore("Computing Cube scores");
    for (auto&& cube : cubeMap_) {
        if (config_->performDiagonalFiltering()) {    // do not filter if M4 parameters are not set
            filterLinksInCube(cube.first);
        }
        if (config_->oldCubeScore()) {
            scoreToCube_[computeCubeScoreOld(cube.first)].emplace(cube.first);
        } else {
            scoreToCube_[computeCubeScore(cube.first)].emplace(cube.first);
        }
    }
    tsScore.endAndPrint(Timestep::minutes);
    std::cout << "Memory usage after Score computation" << std::endl << mm << std::endl;
}



size_t Cubeset::area(Cube const & cube) const {
    auto idMap = linkset_->idMapping();
    auto s0id = cube.tiledistance(0).sequence();
    auto s0len = static_cast<long long>(fastaCollection_->sequenceLength(s0id, *idMap));
    auto s1id = cube.tiledistance(1).sequence();
    auto s1len = static_cast<long long>(fastaCollection_->sequenceLength(s1id, *idMap));
    auto dist = cube.tiledistance(1).distance();
    auto F = static_cast<long long>(config_->tileSize());

    auto v = std::min(s1len - (F*dist), s0len);
    auto u = std::max(-1*F*dist, 0ll);
    auto A = (v - u) * F;
    return std::max(A, 0ll);    // estimate may be negative
}



double Cubeset::computeCubeScore(std::shared_ptr<Cube const> cube) const {
    if (cube->dimensionality() != 2) { std::cerr << "[WARNING] -- Cubeset::computeCubeScore -- CURRENTLY ONLY SUPPORT EXACTLY TWO INPUT GENOMES, SCORE WILL BE WRONG! (Currently " << cube->dimensionality() << " genomes in cube)" << std::endl; }
    auto subtileFrequency = std::vector<size_t>(cubeScoreParameter_, 0); // create each bucket with count 0
    auto& linksInCube = cubeMap_.at(cube);
    for (auto&& link : linksInCube) {
        auto tileindex = cube->tileindex(1);
        auto i = static_cast<long long>(link->position(0));
        auto j = static_cast<long long>(link->position(1));
        auto tilestart = tileindex * static_cast<long long>(config_->tileSize());
        auto jscaled = static_cast<double>(j - i - tilestart);
        auto bucket = static_cast<size_t>( std::floor( jscaled * static_cast<double>(cubeScoreParameter_) / static_cast<double>(config_->tileSize()) ) );
        if (jscaled < 0) { throw std::runtime_error("[ERROR] -- Cubeset::computeCubeScore -- jscaled is negative ((i,j) = (" + std::to_string(link->position(0)) + ", " + std::to_string(link->position(1)) + ") resp. (" + std::to_string(i) + ", " + std::to_string(j) + "))");}
        subtileFrequency.at(bucket) += 1;//linkCount;
    }
    double score = 0;
    auto A = static_cast<double>(std::max(area(*cube), cubeAreaCutoff_));
    for (auto&& h : subtileFrequency) {
        auto n = static_cast<double>(h*h);
        n *= static_cast<double>(eta_)/A;
        score += n;
    }
    if (score < 0) { throw std::runtime_error("[ERROR] -- Cubeset::computeCubeScore -- negative score"); }
    return score;
}



double Cubeset::computeCubeScoreOld(std::shared_ptr<Cube const> cube) const {
    auto linkDelta = [](std::shared_ptr<Link const> const & lhs, std::shared_ptr<Link const> const & rhs) {
        double delta;
        if (lhs->dimensionality() == 2) {
            auto diagLhs = static_cast<long long>(lhs->occurrence(0).position())
                             - static_cast<long long>(lhs->occurrence(1).position());
            auto diagRhs = static_cast<long long>(rhs->occurrence(0).position())
                             - static_cast<long long>(rhs->occurrence(1).position());
            delta = static_cast<double>(std::abs(diagRhs - diagLhs));
        } else {
            std::vector<size_t> lhsStartPoint(lhs->dimensionality());
            std::vector<size_t> rhsStartPoint(rhs->dimensionality());
            size_t lhsMin = ULLONG_MAX;
            size_t rhsMin = ULLONG_MAX;
            for (size_t i = 0; i < lhs->dimensionality(); ++i) {
                lhsStartPoint.at(i) = lhs->occurrence(i).position();
                rhsStartPoint.at(i) = rhs->occurrence(i).position();
                lhsMin = (lhsStartPoint.at(i) < lhsMin) ? lhsStartPoint.at(i) : lhsMin;
                rhsMin = (rhsStartPoint.at(i) < rhsMin) ? rhsStartPoint.at(i) : rhsMin;
            }
            size_t squareSum = 0;
            for (size_t i = 0; i < lhs->dimensionality(); ++i) {
                auto lhsPos = lhsStartPoint.at(i) - lhsMin;
                auto rhsPos = rhsStartPoint.at(i) - rhsMin;
                squareSum += (lhsPos < rhsPos)
                        ? std::pow((rhsPos - lhsPos), 2)
                        : std::pow((lhsPos - rhsPos), 2);
            }
            delta = std::sqrt(squareSum);
        }
        return delta;
    };

    double mu = static_cast<double>(cubeScoreParameter_);
    auto& linksInCube = cubeMap_.at(cube);
    // single-Link Cubes
    if (linksInCube.size() == 1) {
        auto linkCount = linkset_->linkCount(*linksInCube.begin());
        if (linkCount == 0) { throw std::runtime_error("[ERROR] -- Cubeset::computeCubeScoreOld -- Linkcount is zero"); }
        if (linkCount == 1) {
            return 1.;
        } else {
            return static_cast<double>(linkCount) * (1. + mu);
        }
    }
    if (!std::numeric_limits<double>::has_infinity) {
        throw std::runtime_error("[ERROR] -- Cubeset::computeCubeScoreOld -- Type `double` cannot represent infinity on this machine");
    }
    auto links = std::vector<std::shared_ptr<Link const>>{linksInCube.begin(), linksInCube.end()};
    std::sort(links.begin(), links.end(), LinkPtrLess());
    double score = 0;
    auto previousDelta = std::numeric_limits<double>::infinity();
    auto it = links.begin();
    auto nextIt = std::next(links.begin(), 1);
    while (it != links.end()) {
        auto delta = (nextIt == links.end())
                ? std::numeric_limits<double>::infinity()
                : linkDelta(*it, *nextIt);
        auto linkCount = linkset_->linkCount(*it);
        double thisScore = 0;
        if (linkCount > 1) {
            // if a Link was created multiple times, act as if it was there this number of times
            //   i.e. add the multiplied scores for the in-between identical diagonals
            //   as the algorithm would definitely choose the same-diagonal links for scoring
            thisScore = static_cast<double>(linkCount) * (mu + 1.);
        } else {
            // otherwise, compare the neighbouring deltas and compute score properly
            auto minDelta = std::min(previousDelta, delta);
            thisScore = 1. + (mu / (1. + minDelta));
        }
        score += thisScore;
        previousDelta = delta;
        ++it;
        ++nextIt;
    }

    // normalization
    auto A = static_cast<double>(std::max(area(*cube), cubeAreaCutoff_));
    score *= static_cast<double>(eta_)/A;

    return score;
}




void Cubeset::filterLinksInCube(std::shared_ptr<Cube const> cube) {
    auto & linksInCube = cubeMap_.at(cube);
    std::vector<Link> filterableLinks;
    for (auto&& linkptr : linksInCube) { filterableLinks.emplace_back(*linkptr); }
    linksInCube.clear(); // delete Links from set
    std::sort(filterableLinks.begin(), filterableLinks.end());
    std::vector<Link> result;
    auto filter = DiagonalMatchesFilter<Link>(config_);
    filter.applyDiagonalMatchesFilter(filterableLinks, result, true);   // true: quiet (do not print Timesteps)
    for (auto&& link : result) { linksInCube.insert(std::make_shared<Link>(link)); }
}
