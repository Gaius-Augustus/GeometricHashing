#include "Cubeset.h"

Cubeset::Cubeset(std::shared_ptr<Linkset const> linkset,
                 std::shared_ptr<Configuration const> config)
    : config_{config}, cubeMap_{},
      linkset_{linkset}, mu_{config->cubeScoreMu()},
      numGenomes_{linkset->numGenomes()}, scoreToCube_{}/*, tiledistToCube_{}*/ {
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

    Timestep tsScore("Computing only Cube scores");
    for (auto&& cube : cubeMap_) {
        if (config_->performDiagonalFiltering()) {    // do not filter if M4 parameters are not set
            filterLinksInCube(cube.first);
        }
        scoreToCube_[computeCubeScore(cube.first)].emplace(cube.first);
    }
    tsScore.endAndPrint(Timestep::minutes);
    std::cout << "Memory usage after Score computation" << std::endl << mm << std::endl;
}



double Cubeset::computeCubeScore(std::shared_ptr<Cube const> cube) const {
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

    auto& linksInCube = cubeMap_.at(cube);

    // single-Link Cubes
    if (linksInCube.size() == 1) {
        auto linkCount = linkset_->linkCount(*linksInCube.begin());
        if (linkCount == 0) { throw std::runtime_error("[ERROR] -- Cubeset::computeCubeScore -- Linkcount is zero"); }
        if (linkCount == 1) {
            return 1.;
        } else {
            return static_cast<double>(linkCount) * (1. + mu_);
        }
    }

    if (!std::numeric_limits<double>::has_infinity) {
        throw std::runtime_error("[ERROR] -- Cubeset::computeCubeScore -- Type `double` cannot represent infinity on this machine");
    }

    double score = 0;
    auto previousDelta = std::numeric_limits<double>::infinity();
    auto it = linksInCube.begin();
    auto nextIt = std::next(linksInCube.begin(), 1);
    while (it != linksInCube.end()) {
        auto delta = (nextIt == linksInCube.end())
                ? std::numeric_limits<double>::infinity()
                : linkDelta(*it, *nextIt);
        auto linkCount = linkset_->linkCount(*it);
        if (linkCount > 1) {
            // if a Link was created multiple times, act as if it was there this number of times
            //   i.e. add the multiplied scores for the in-between identical diagonals
            //   as the algorithm would definitely choose the same-diagonal links for scoring
            score += static_cast<double>(linkCount) * (mu_ + 1.);
        } else {
            // otherwise, compare the neighbouring deltas and compute score properly
            auto minDelta = std::min(previousDelta, delta);
            auto thisScore = 1. + (mu_ / (1. + minDelta));
            score += thisScore;
        }
        previousDelta = delta;
        ++it;
        ++nextIt;
    }
    return score;
}



void Cubeset::filterLinksInCube(std::shared_ptr<Cube const> cube) {
    auto & linksInCube = cubeMap_.at(cube);
    std::set<Link> filterableLinks;
    for (auto&& linkptr : linksInCube) { filterableLinks.emplace(*linkptr); }
    linksInCube.clear(); // delete Links from set
    std::vector<Link> result;
    auto filter = DiagonalMatchesFilter<Link>(config_);
    filter.applyDiagonalMatchesFilter(filterableLinks, result, true);   // true: quiet (do not print Timesteps)
    for (auto&& link : result) { linksInCube.insert(std::make_shared<Link>(link)); }
}
