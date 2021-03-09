#ifndef REGIONTUPLEEXTRACTION_H
#define REGIONTUPLEEXTRACTION_H

#include <math.h>
#include <memory>
#include <vector>

#include <tsl/hopscotch_set.h>
#include "Link.h"

// Expects a link container, does not need to be sorted
//template <typename Linkcontainer>
inline auto regionTupleExtraction(//Linkcontainer const & links,
                                  tsl::hopscotch_set<LinkPtr, LinkPtrHash> const & links,
                                  size_t referenceSequenceLength, size_t tilesize,
                                  double alpha, double beta, double gamma) {
    static_assert(std::numeric_limits<double>::is_iec559, "[ERROR] -- regionTupleExtraction -- Type `double` cannot represent negative infinity on this machine");
    // initialize
    auto ntiles = (tilesize) ? static_cast<size_t>(std::floor(static_cast<double>(referenceSequenceLength)/static_cast<double>(tilesize)) + 1) : 1; // if tilesize==0, everything is in one single tile
    auto A = std::vector<double>(ntiles);
    auto Am1 = std::vector<bool>(ntiles); // track inner max of A
    auto B = std::vector<double>(ntiles);
    auto Aprev = -std::numeric_limits<double>::infinity();
    double Bprev = 0;
    double nprev = 0;
    // count links per tile
    auto nlinks = std::vector<size_t>(ntiles, 0);
    for (auto& link : links) {
        if (link->first().genome() == 0) {
            auto tile = (tilesize) ? static_cast<size_t>(std::floor(static_cast<double>(link->position(0))/static_cast<double>(tilesize))) : 0; // if tilesize==0, everythin is in zero-th tile
            try {
                ++(nlinks.at(tile));
            } catch (std::exception const & e) {
                std::cerr << "[ERROR] -- regionTupleExtraction() -- increment link count for tile " << tile << " (with " << ntiles << " ntiles and link position " << link->position(0) << ")" << std::endl;
                throw e;
            }
        }
    }
    // forward
    for (size_t n = 0; n < ntiles; ++n) {
        auto innerA1 = Bprev - beta - gamma;
        auto innerA2 = Aprev - (beta*(static_cast<double>(n)-nprev));
        Am1.at(n) = (innerA1 > innerA2);
        A.at(n) = alpha * static_cast<double>(nlinks.at(n)) + std::max(innerA1, innerA2);
        B.at(n) = std::max(A.at(n), Bprev);
        Aprev = A.at(n);
        Bprev = B.at(n);
        nprev = (nlinks.at(n)) ? n : nprev;
    }
    // backward
    std::vector<std::array<size_t, 2>> intervals;
    bool intervalOpen = false;
    for (int j = ntiles-1; j >= 0; --j) {
        // choose A or B
        if (A.at(j) >= B.at(j)) {
            // A
            if (!intervalOpen) {
                intervals.emplace_back(std::array<size_t, 2>{0,static_cast<size_t>(j)});
                intervalOpen = true;
            }
            if (Am1.at(j)) {
                intervals.rbegin()->at(0) = j;
                intervalOpen = false;
            }
        } else if (intervalOpen) {
            // B and open interval
            intervals.rbegin()->at(0) = j;
            intervalOpen = false;
        }
    }
    if (intervalOpen) {
        intervals.rbegin()->at(0) = 0;
    }
    return intervals;
}

#endif // REGIONTUPLEEXTRACTION_H
