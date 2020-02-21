#ifndef LINK_H
#define LINK_H

#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

#include "prettyprint/prettyprint.hpp"
#include "CustomHashGeneral.h"
#include "KmerOccurrence.h"

//! Representation of a link
/*! Connects a set of tiles across several genomes,
 * only one tile per genome is allowed and the
 * reference genome must be present */
class Link {
public:
    //! Constructor, initializes empty object
    Link() : occurrence_{} {}
    //! Constructor (2), create Link from vector of KmerOccurrence s
    /*! \param occurrences Vector of KmerOccurrence s
     *
     * \details Sorts the \c tile_ vector if it is not sorted yet,
     * throws if more than one occurrence per genome is in \c occurrences */
    Link(std::vector<KmerOccurrence> const & occurrences)
        : occurrence_{occurrences.begin(), occurrences.end()} {
        // check sorted as well as each genome at most once
        auto first = occurrence_.begin();
        auto last = occurrence_.end();
        if (first != last) {
            std::unordered_set<size_t> seenGenome;
            bool sorted = true;
            seenGenome.emplace(first->genome());
            auto next = first;
            while (++next != last) {    // http://www.cplusplus.com/reference/algorithm/is_sorted/
                if (seenGenome.find(next->genome()) != seenGenome.end()) {
                    throw std::runtime_error("[ERROR] -- Link::Link (2) -- Only one occurrence per genome allowed");
                }
                seenGenome.emplace(next->genome());
                if (*next < *first) {
                    sorted = false;
                }
                ++first;
            }
            if (!sorted) {
                std::sort(occurrence_.begin(), occurrence_.end());
            }
        }
    }

    //! Returns the number of genomes that appear in this Link
    auto dimensionality() const { return static_cast<uint16_t>(occurrence_.size()); }
    //! Getter for the i-th genome ID in this Link
    auto genome(uint32_t i) const { return occurrence_.at(i).genome(); }
    //! Append a tile to this Link, expensive as this involves sorted insert and check if only one occurrence per genome
    void insertOccurrence(uint8_t genomeID,  uint32_t sequenceID, size_t tileID, bool reverse, std::string const & kmer) {
        insertOccurrence(KmerOccurrence(genomeID, sequenceID, tileID, reverse, kmer));
    }
    void insertOccurrence(KmerOccurrence const & newTile) {
        if (occurrence_.size() >= UINT32_MAX) {
            auto message = "[ERROR] -- Link::insertTile -- Link cannot hold more than " + std::to_string(UINT16_MAX) + " tiles";
            throw std::runtime_error(message);
        }
        for (auto&& occ : occurrence_) {
            if (occ.genome() == newTile.genome()) { throw std::runtime_error("[ERROR] -- Link::insertTile -- Only one occurrence per genome allowed"); }
        }
        occurrence_.insert(std::upper_bound(occurrence_.begin(), occurrence_.end(), newTile),
                     newTile);
    }
    //! Getter for the internal occurrence vector
    auto const & occurrence() const { return occurrence_; }
    //! Get i-th KmerOccurrence in this Link
    auto const & occurrence(size_t i) const { return occurrence_.at(i); }
    //! Implements operator== for Link s by checking if tile vectors and kmer strings are equal
    bool operator==(Link const & rhs) const {
        return std::equal(occurrence_.begin(), occurrence_.end(),
                          rhs.occurrence_.begin(), rhs.occurrence_.end(),
                          KmerOccurrence::equalSpot);
    }
    //! Compares two Link s if one is 'less' than the other
    /*! Less is determined in the following manner:
     * First of all, Occurrences in Link s are sorted and as there is only one
     * occurrence per genome allowed, they are strictly sorted after the genome IDs.
     *
     * In a first step, a decision is sought using the occurrences in the shared genomes,
     * iterating through the shared genomes. As soon as the Link s differ there, we're done.
     *
     * If all shared genome occurrences are at the same spot, a decision is sought based
     * on dimensionality, i.e. the Link with fewer Occurrences is less.
     *
     * If the Link s have equal dimensionality, a decision is sought by looking at which
     * Link has the lowest genome occurrrence in the non-shared genomes
     *
     * If a decision could not be made until here, it means that they have the same occurrences
     * in the shared genomes, have the same number of genomes and there are no not-shared
     * genomes, i.e. the Link s are equal, thus this operator returns false. */
    bool operator<(Link const & rhs) const {
        auto thisTileIt = occurrence_.begin();
        auto rhsTileIt = rhs.occurrence_.begin();
        while(thisTileIt != occurrence_.end() && rhsTileIt != rhs.occurrence_.end()) {
            if (thisTileIt->genome() == rhsTileIt->genome()) {  // if Links share genome, check occurrence
                if ((*thisTileIt) < (*rhsTileIt)) {
                    return true;    // less
                } else if (!KmerOccurrence::equalSpot(*thisTileIt, *rhsTileIt)) {
                    return false;   // greater
                }
            }
            // if occurrences are equal, find next shared genome
            ++thisTileIt;   // next tile in this
            if (thisTileIt == occurrence_.end()) { break; }   // UB otherwise
            while (thisTileIt->genome() > rhsTileIt->genome()) {
                ++rhsTileIt;    // while rhs genome id is lower than this genome id, iterate rhs
                if (rhsTileIt == rhs.occurrence_.end()) { break; }   // UB otherwise
            }
        }
        // no shared occurrences differ, check dimensionality
        if (occurrence_.size() < rhs.occurrence_.size()) {
            return true;
        } else if (occurrence_.size() > rhs.occurrence_.size()) {
            return false;
        }
        // no shared occurrences differ and equal dimensionality, check lowest not-shared genome id
        for (size_t i = 0; i < occurrence_.size(); ++i) {
            if (occurrence_.at(i).genome() < rhs.occurrence_.at(i).genome()) {
                return true;
            } else if (occurrence_.at(i).genome() > rhs.occurrence_.at(i).genome()) {
                return false;
            }
        }
        return false;
    }
    //! Implements operator<< for a Link object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, Link const & l) {
        out << l.occurrence_ << std::endl;
        return out;
    }
    //! Getter for the i-th position in this Link
    auto position(uint32_t i) const { return occurrence_.at(i).position(); }
    //! Getter for the i-th strand information in this Link
    auto reverse(uint32_t i) const { return occurrence_.at(i).reverse(); }
    //! Getter for the i-th sequence ID in this Link
    auto sequence(uint32_t i) const { return occurrence_.at(i).sequence(); }

    //! Returns the first occurrence in this Link, throws if Link is empty
    auto const & first() const { return occurrence(0); }
    //! Returns the second occurrence in this Link, throws if Link has less than two occurrences
    auto const & second() const { return occurrence(1); }
    //! True if the first two occurrences of this and rhs have the same distance, respectively
    auto sameDistance(Link const & rhs) const {
        auto pair = KmerOccurrencePair(first(), second());
        auto rhsPair = KmerOccurrencePair(rhs.first(), rhs.second());
        return pair.sameDistance(rhsPair);
    }

private:
    //! Stores the tiles in this link
    std::vector<KmerOccurrence> occurrence_;
};



//! Computes hash value of a Link object
struct LinkPtrHash {
    size_t operator()(std::shared_ptr<Link const> const & link) const {
        KmerOccurrencePositionHash hashfun;
        size_t seed = 0;
        for (auto&& occ : link->occurrence()) {
            customCombineHash(seed, hashfun(occ));
        }
        return seed;
    }
};



//! Implements operator() to check if two \c std::shared_ptr<Link> point to equal Link s
struct LinkPtrEqual {
    bool operator()(std::shared_ptr<Link const> const & lhs, std::shared_ptr<Link const> const & rhs) const {
        return *lhs == *rhs;
    }
};



//! Implements operator() to check if first \c std::shared_ptr<Link> is less than the second w.r.t. the pointed-to Link s
struct LinkPtrLess {
    bool operator()(std::shared_ptr<Link const> const & lhs, std::shared_ptr<Link const> const & rhs) const {
        return *lhs < *rhs;
    }
};



//! Overloads operator== for std::shared_ptr<Link const> for easier set comparison
inline bool operator==(std::shared_ptr<Link const> const & lhs, std::shared_ptr<Link const> const & rhs) {
    return *lhs == *rhs;
}



//! Implements operator<< for a \c std::shared_ptr of Link for more convenient \c std::ostream functionality
inline std::ostream& operator<<(std::ostream& out, std::shared_ptr<Link const> const & lp) {
    out << *lp;
    return out;
}

#endif // LINK_H
