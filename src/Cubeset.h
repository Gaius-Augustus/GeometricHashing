#ifndef CUBESET_H
#define CUBESET_H

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <fstream>  // for writing statistics to file
#include <iostream>
#include <iterator>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>
#include "json/json.hpp"
#include "prettyprint/prettyprint.hpp"
#include "Cube.h"
#include "DiagonalMatchesFilter.h"
#include "Linkset.h"
#include "MemoryMonitor.h"
#include "ProgressBar.h"
#include "Timestep.h"



//! Sort Link s from the same Cube after their diagonals (like in KmerOccurrencePair)
/*! Undefined behaviour if the Link s are from different Cube s, i.e. if
 * genomes, sequences and strands are not identical in each respective KmerOccurrence */
class LinkInCubeDiagonalCompare {
public:
    bool operator()(std::shared_ptr<Link const> const & lhs, std::shared_ptr<Link const> const & rhs) const {
        for (uint16_t i = 1; i < lhs->dimensionality(); ++i) {
            // check genomes in increasing ID order and decide on the absolute position distance of the occurrences in the respective genome
            auto lhsDistance = static_cast<long long>(lhs->occurrence(0).position()) - static_cast<long long>(lhs->occurrence(i).position());
            auto rhsDistance = static_cast<long long>(rhs->occurrence(0).position()) - static_cast<long long>(rhs->occurrence(i).position());
            if (lhsDistance < rhsDistance) {
                return true;
            } else if (lhsDistance > rhsDistance) {
                return false;
            }
        }
        // same diagonal at this point, sort after absolute positions
        for (uint16_t i = 0; i < lhs->dimensionality(); ++i) {
            if (lhs->occurrence(i).position() < rhs->occurrence(i).position()) {
                return true;
            } else if (lhs->occurrence(i).position() > rhs->occurrence(i).position()) {
                return false;
            }
        }
        return *lhs < *rhs; // do normal comparison
    }
};



//! Represents the set of all Cube s
/*! Cubes are generated from a Linkset and a mapping from
 * Cubes to a set of Link s that is contained in the respective cube,
 * as well as a mapping from Tiledistance s to the cubes that share these
 * distances is stored. Finally, artificial cubes and the relation between
 * cubes (Hasse diagram) are computed and stored. */
class Cubeset {
public:
    //! Factory function that computes the score of each cube
    /*! Executed during Hasse computation but after all respective
     * predecessors are known */
    double computeCubeScore(std::shared_ptr<Cube const> cube) const;

    //! Constructor
    /*! \param linkset Linkset from which the Cubeset is computed
     *
     * \details Iterates over all Links s and creates the respective Cube s,
     * stores the mappings from Cube to set of Link s and from Tiledistance
     * to set of Cube s */
    Cubeset(std::shared_ptr<Linkset const> linkset,
            std::shared_ptr<Configuration const> config);

    //! Getter for member \c cubeMap_
    auto const & cubeMap() const { return cubeMap_; }
    //! Returns \c std::unordered_set of Link s in \c cube
    auto const & links(std::shared_ptr<Cube const> const & cube) const { return cubeMap_.at(cube); }
    //! Getter for member \c scoreToCube_
    auto const & scoreToCube() const { return scoreToCube_; }
    //! Provide const access to the underlying Linkset
    auto const & underlyingLinkset() const { return linkset_; }
    //! Implements operator<< for a Cubeset object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, Cubeset const & cs) {
        out << "Cube => Set of Links" << std::endl;
        for (auto&& element : cs.cubeMap_) {
            out << element.first << " => " << element.second << std::endl;
        }
        out << "Score => Cubes" << std::endl;
        for (auto&& element : cs.scoreToCube_) {
            out << element.first << " => " << element.second << std::endl;
        }
        return out;
    }

private:
    //! M4 filtering
    void filterLinksInCube(std::shared_ptr<Cube const> cube);

    //! Configuration for filtering
    std::shared_ptr<Configuration const> config_;
    //! Stores the mapping from each Cube to a tuple of set of links and vector of direct predecessors
    /*! Implicitly stores the set of all Cube s as keys */
    std::map<std::shared_ptr<Cube const>,
             std::set<std::shared_ptr<Link const>, // Links in Cube
                      LinkInCubeDiagonalCompare>,
             CubePtrLess> cubeMap_;
    //! Linkset from which this Cubeset was created
    std::shared_ptr<Linkset const> const linkset_;
    //! Cube scoring parameter mu
    double mu_;
    //! Number of genomes in input data
    size_t numGenomes_;
    //! Stores a score to Cube (s) mapping for all Cube s
    std::map<double, std::unordered_set<std::shared_ptr<Cube const>,
                                        CubePtrHash,
                                        CubePtrEqual>> scoreToCube_;
};

#endif // CUBESET_H
