#ifndef CUBE_H
#define CUBE_H

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>
#include "prettyprint/prettyprint.hpp"
#include "Link.h"
#include "MemoryMonitor.h"



//! Representation of a Distance between two genome tiles
/*! Basically the same functionality as a KmerOccurrence but since
 * that can't store negative 'positions' (in this case, tile distances),
 * this class takes care of the neccessary modifications */
class Tiledistance {
public:
    //! Constructor
    Tiledistance(uint8_t genomeID, uint32_t sequenceID, long long distance, bool reverseStrand)
        : occurrence_{KmerOccurrence(0,0,0,false,"")} {
        // bits for position: pppp pppp pppp pppp pppp pppp pppp pppp pppp pppp
        // i.e.             0xf    f    f    f    f    f    f    f    f    f
        // LSB needed for sign, so only allow numbers [-0x7fffffffff, 0x7fffffffff]
        auto negative = false;
        if (distance < 0) {
            negative = true;
            distance *= -1;
        }
        if (distance > 0x7fffffffff) { throw std::runtime_error("[ERROR] -- Tiledistance -- Distance too large"); }
        auto bits = std::bitset<40>(distance);
        bits <<= 1;
        if (negative) { bits.set(0); }
        occurrence_ = KmerOccurrence(genomeID, sequenceID, bits.to_ullong(), reverseStrand, "");
    }
    auto distance() const {
        auto bits = std::bitset<80>(occurrence_.position());
        auto negative = bits.test(0);
        bits >>= 1; // remove sign bit
        auto totalDist = bits.to_ullong();
        auto distance = (negative)
                ? static_cast<long long>(totalDist) * -1
                : static_cast<long long>(totalDist);
        return distance;
    }
    auto genome() const { return occurrence_.genome(); }
    auto const & kmerOccurrence() const { return occurrence_; }
    bool operator<(Tiledistance const & rhs) const { return occurrence_ < rhs.occurrence_; }
    bool operator==(Tiledistance const & rhs) const { return occurrence_ == rhs.occurrence_; }
    auto reverse() const { return occurrence_.reverse(); }
    auto sequence() const { return occurrence_.sequence(); }
    friend std::ostream& operator<<(std::ostream& out, Tiledistance const & td) {
        //out << td.occurrence_; // prints the wrong distances
        out << "(" << static_cast<uint>(td.genome()) << ", " << td.sequence() << ", "
            << td.distance() << ", " << td.reverse() << ")";
        return out;
    }
private:
    KmerOccurrence occurrence_;
};



//! Representation of a cube
/*! A cube is basically a \c vector of KmerOccurrence s that resemble tiles relative to a reference position.*/
class Cube {
public:
    //! Constructor (1)
    /*! \param link Link from which the Cube is created
     * \param tileSize Size of tiles
     *
     * \details Computes and stores the Tiledistance of each occurrence (including reference)
     * in the Link w.r.t. the reference occurrence. Two identical Cube s can be created
     * from different Link s. The Tiledistance s have the same order as the
     * occurrences in \c link. */
    Cube(std::shared_ptr<Link const> link, size_t tileSize)
        : tiledistance_{} {
        if (link->genome(0) != 0) { throw std::runtime_error("[ERROR] -- Cube::Cube -- Cannot create Cubes from Links that miss reference genome (0)"); }
        auto refPosition = link->position(0);
        for (auto&& occurrence : link->occurrence()) {
            tiledistance_.emplace_back(occurrence.genome(),
                                       occurrence.sequence(),
                                       positionsToTile(occurrence.position(), refPosition, tileSize),
                                       occurrence.reverse());
        }
    }

    //! Returns the number of genomes present in this Cube
    uint16_t dimensionality() const { return static_cast<uint16_t>(tiledistance_.size()); }
    //! Getter for the tiledistance_ member
    auto const & tiledistance() const {
        return tiledistance_;
    }
    //! Getter for the i-th Tiledistance
    auto const & tiledistance(uint16_t i) const {
        return tiledistance_.at(i);
    }
    //! Getter for i-th tile index
    auto tileindex(uint16_t i) const {
        return tiledistance_.at(i).distance();
    }
    //! Operator== for Cube
    bool operator==(Cube const & rhs) const {
        return tiledistance_ == rhs.tiledistance_;
    }
    //! Operator< for Cube
    bool operator<(Cube const & rhs) const {
        return tiledistance_ < rhs.tiledistance_;
    }
    //! Computes floor((j-i)/F) where i is position is reference genome and F is tilesize
    static long long positionsToTile(size_t position, size_t referencePosition, size_t tileSize) {
        auto j = static_cast<double>(position);
        auto i = static_cast<double>(referencePosition);
        return static_cast<long long>( std::floor((j-i)/static_cast<double>(tileSize)) );
    }
    friend std::ostream& operator<<(std::ostream& out, Cube const & c) {
        out << c.tiledistance_;
        return out;
    }

private:
    //! Stores the Tiledistance s
    std::vector<Tiledistance> tiledistance_;
};



//! Computes hash value of a Cube object
struct CubePtrHash {
    size_t operator()(std::shared_ptr<Cube const> const & cp) const {
        KmerOccurrencePositionHash hashfun;
        size_t seed = 0;
        for (auto&& occ : cp->tiledistance()) {
            customCombineHash(seed, hashfun(occ.kmerOccurrence()));
        }
        return seed;
    }
};



//! Implements operator() to check if two \c std::shared_ptr<Cube> point to equal Cube s
struct CubePtrEqual {
    bool operator()(std::shared_ptr<Cube const> const & lhs,std::shared_ptr<Cube const> const & rhs) const {
        return *lhs == *rhs;
    }
};



//! Implements operator() to check if first \c std::shared_ptr<Cube> point to a 'less' Cube than the second
struct CubePtrLess {
    bool operator()(std::shared_ptr<Cube const> const & lhs, std::shared_ptr<Cube const> const & rhs) const {
        return *lhs < *rhs;
    }
};



//! Implements operator<< for a \c std::shared_ptr of Cube for more convenient \c std::ostream functionality
inline std::ostream& operator<<(std::ostream& out, std::shared_ptr<Cube const> const & cp) {
    out << *cp;
    return out;
}

#endif // CUBE_H
