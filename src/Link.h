#ifndef LINK_H
#define LINK_H

#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

#include "prettyprint.hpp"
#include "CustomHashGeneral.h"
#include "KmerOccurrence.h"



//! Reusable template function to combine overlapping matches in a sorted container
template <typename MatchContainer>
void groupOverlappingSeeds(MatchContainer & sortedMatches) {
    // copy matches, clear input container
    std::vector<typename MatchContainer::value_type> matches(sortedMatches.begin(), sortedMatches.end());
    sortedMatches.clear();
    // combine overlapping seeds
    auto lIt = matches.begin();
    if (lIt != matches.end()) {
        auto rIt = std::next(matches.begin()); // begin() + 1
        while (lIt != matches.end()) {
            auto rightBorder = lIt->first().position() + lIt->span() - 1;
            size_t extend = 0;
            while (rIt != matches.end() && lIt->sameDiagonal(*rIt) && rIt->first().position() <= rightBorder) {
                auto newRightBorder = rIt->first().position() + rIt->span() - 1;
                if (newRightBorder > rightBorder) {
                    extend += newRightBorder - rightBorder;
                    rightBorder = newRightBorder;
                }
                ++rIt;
            }
            lIt->extendSpanToRight(extend);
            sortedMatches.insert(sortedMatches.end(), *lIt);    // re-add extended match
            lIt = rIt; // next match group
            if (rIt != matches.end()) { ++rIt; }
        }
    }
}



//! compute chunk ID from position sum
inline size_t chunkID_impl(size_t sum, size_t chunksize) {
    return static_cast<size_t>(std::floor(static_cast<long double>(sum) / static_cast<long double>(chunksize)));
}



//! Representation of a link
/*! Connects a set of tiles across several genomes,
 * only one tile per genome is allowed and the
 * reference genome must be present */
class Link {
public:
    //! Constructor, initializes empty object
    Link() : occurrence_{}, span_{1} {}
    //! Constructor (2), create Link from vector of KmerOccurrence s
    /*! \param occurrences Vector of KmerOccurrence s
     *
     * \details Sorts the \c tile_ vector if it is not sorted yet,
     * throws if more than one occurrence per genome is in \c occurrences */
    Link(std::vector<KmerOccurrence> const & occurrences, size_t span)
        : occurrence_{occurrences.begin(), occurrences.end()}, span_{span} {
        if (span_ < 1) { throw std::runtime_error("[ERROR] -- Link -- span must be at least 1"); }
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

    //! Return chunk ID of \c this
    size_t chunkID(size_t chunksize) const {
        if (chunksize == 0) { return 0; }
        size_t sum = 0;
        for (auto&& occ : occurrence_) { sum += occ.position(); }
        return chunkID_impl(sum, chunksize);
    }
    //! Returns the diagonal of this link (i.e. a vector [i-i, j-i, k-i, ...] for each occurrence i,j,k,...)
    auto diagonal() const {
        std::vector<long long> diag;
        for (size_t i = 0; i < occurrence_.size(); ++i) {
            diag.emplace_back(static_cast<long long>(occurrence_[i].position())
                              - static_cast<long long>(occurrence_[0].position()));
        }
        return diag;
    }
    //! Returns the number of genomes that appear in this Link
    auto dimensionality() const { return occurrence_.size(); }
    //! Possible to extend span to the right by \c amount (i.e. span += amount)
    /*! e.g. if there is a neighbouring Link on the same diagonal
     *   that should be combined with this Link */
    void extendSpanToRight(size_t amount) { span_ += amount; }
    //! Getter for the i-th genome ID in this Link
    auto genome(uint32_t i) const { return occurrence_.at(i).genome(); }
    //! Append a tile to this Link, expensive as this involves sorted insert and check if only one occurrence per genome
    void insertOccurrence(uint8_t genomeID,  uint32_t sequenceID, size_t tileID, bool reverse, std::string const & kmer) {
        insertOccurrence(KmerOccurrence(genomeID, sequenceID, tileID, reverse, kmer));
    }
    void insertOccurrence(KmerOccurrence const & newTile) {
        for (auto&& occ : occurrence_) {
            if (occ.genome() == newTile.genome()) { throw std::runtime_error("[ERROR] -- Link::insertOccurrence -- Only one occurrence per genome allowed"); }
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
                          KmerOccurrence::equalSpot)
                && span_ == rhs.span_;
    }
    //! Compares two Link s if one is 'less' than the other
    bool operator<(Link const & rhs) const {
        // first level: compare link dimensions
        if (occurrence_.size() == rhs.occurrence_.size()) {
            if (occurrence_.size() == 0) { return false; } // empty links
            // check each genome, as soon as genome/sequence/strand differ, return
            for (size_t i = 0; i < occurrence_.size(); ++i) {
                if (occurrence_[i].genome() != rhs.occurrence_[i].genome()) {
                    return occurrence_[i].genome() < rhs.occurrence_[i].genome();
                }
                if (occurrence_[i].sequence() != rhs.occurrence_[i].sequence()) {
                    return occurrence_[i].sequence() < rhs.occurrence_[i].sequence();
                }
                if (occurrence_[i].reverse() != rhs.occurrence_[i].reverse()) {
                    return occurrence_[i].reverse() < rhs.occurrence_[i].reverse();
                }
            }
            // at this point, links are from same sequences (and strands) -> sort diagonals together
            if (diagonal() != rhs.diagonal()) {
                return diagonal() < rhs.diagonal();
            }
            // same diagonal, sufficient to compare first position
            if (occurrence_[0].position() != rhs.occurrence_[0].position()) {
                return occurrence_[0].position() < rhs.occurrence_[0].position();
            }
            // last possibility: differing spans
            return span_ < rhs.span_;
        } else {
            return occurrence_.size() < rhs.occurrence_.size();
        }
    }
    //! Implements operator<< for a Link object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, Link const & l) {
        out << l.occurrence_ << "(span " << l.span_ << ")" << std::endl;
        return out;
    }
    //! Getter for the i-th position in this Link
    auto position(size_t i) const { return occurrence_.at(i).position(); }
    //! Getter for the i-th strand information in this Link
    auto reverse(size_t i) const { return occurrence_.at(i).reverse(); }
    //! True if the same occurrences w.r.t. genome and sequence and same diagonal
    /*! Same diagonal means that the pairwise distances of all occurrences
     *   with \c occurrence_[0] are equal for both \c this and \c rhs.
     *  This implies that all sequences (and thus genomes) in both links are equal */
    bool sameDiagonal(Link const & rhs) const {
        if (occurrence_.size() == rhs.occurrence_.size()) {
            if (occurrence_.size() == 0) { return true; }
            if (!sameSeq(occurrence_[0], rhs.occurrence_[0])) { return false; }
            for (size_t i = 1; i < occurrence_.size(); ++i) {
                if (!sameSeq(occurrence_[i], rhs.occurrence_[i])) { return false; }
                if ((occurrence_[i].position() - occurrence_[0].position())
                        != (rhs.occurrence_[i].position() - rhs.occurrence_[0].position())) {
                    return false;
                }
            }
            return true;
        } else {
            return false;
        }
    }
    //! \c true if both Link s are from the same sequence set
    bool sameSequences(Link const & rhs) const {
        if (occurrence_.size() == rhs.occurrence_.size()) {
            for (size_t i = 0; i < occurrence_.size(); ++i) {
                if (!sameSeq(occurrence_[i], rhs.occurrence_[i])) { return false; }
            }
            return true;
        } else {
            return false;
        }
    }
    //! Getter for the i-th sequence ID in this Link
    auto sequence(size_t i) const { return occurrence_.at(i).sequence(); }
    //! Getter for this Link s span
    auto span() const { return span_; }
    //! Returns the first occurrence in this Link, throws if Link is empty
    auto const & first() const { return occurrence(0); }
    //! Returns the second occurrence in this Link, throws if Link has less than two occurrences
    auto const & second() const { return occurrence(1); }

private:
    //! Helpter function
    bool sameSeq(KmerOccurrence const & lhs, KmerOccurrence const & rhs) const {
        return lhs.genome() == rhs.genome()
                && lhs.sequence() == rhs.sequence()
                && lhs.reverse() == rhs.reverse();
    }
    //! Stores the tiles in this link
    std::vector<KmerOccurrence> occurrence_;
    //! Span (width) of the link
    size_t span_;
};



//! Computes hash of a Link
inline size_t linkhash_impl(Link const & link, size_t seed) {
    KmerOccurrencePositionHash hashfun;
    for (auto&& occ : link.occurrence()) {
        customCombineHash(seed, hashfun(occ));
    }
    return seed;
}
//! Computes hash value of a Link object
struct LinkHash {
    size_t operator()(Link const & link) const {
        return linkhash_impl(link, std::hash<size_t>{}(link.span()));
    }
};
//! Computes hash value of a Link object, ignoring span
struct LinkHashIgnoreSpan {
    size_t operator()(Link const & link) const {
        return linkhash_impl(link, 0);
    }
};

//! Implements operator() to check if two Link s are equal, ignoring span
struct LinkEqualIgnoreSpan {
    bool operator()(Link const & lhs, Link const & rhs) const {
        return std::equal(lhs.occurrence().begin(), lhs.occurrence().end(),
                          rhs.occurrence().begin(), rhs.occurrence().end(),
                          KmerOccurrence::equalSpot);
    }
};



//! Wrapper class around \c std::shared_ptr<Link const>
class LinkPtr {
public:
    LinkPtr() : link_{std::make_shared<Link>()} {}
    LinkPtr(std::vector<KmerOccurrence> const & occurrences, size_t span)
        : link_{std::make_shared<Link>(occurrences, span)} {}
    LinkPtr(Link const & link) : link_{std::make_shared<Link>(link)} {}
    //template<typename... Args>
    //LinkPtr(Args... args) : link_{std::make_shared<Link>(args...)} {}

    size_t chunkID(size_t chunksize) const { return link_->chunkID(chunksize); }
    auto diagonal() const { return link_->diagonal(); }
    auto dimensionality() const { return link_->dimensionality(); }
    auto extendSpanToRight(size_t amount) { link_->extendSpanToRight(amount); }
    auto genome(uint32_t i) const { return link_->genome(i); }
    void insertOccurrence(uint8_t genomeID,  uint32_t sequenceID, size_t tileID, bool reverse, std::string const & kmer) {
        link_->insertOccurrence(genomeID,  sequenceID, tileID, reverse, kmer);
    }
    void insertOccurrence(KmerOccurrence const & newTile) { link_->insertOccurrence(newTile); }
    auto const & link() const { return *link_; }
    auto const & occurrence() const { return link_->occurrence(); }
    auto const & occurrence(size_t i) const { return link_->occurrence(i); }
    bool operator==(LinkPtr const & rhs) const { return *link_ == *(rhs.link_); }
    bool operator<(LinkPtr const & rhs) const { return *link_ < *(rhs.link_); }
    //! direct dereferentiation of underlying ptr
    auto const & operator*() const { return *link_; }
    //! direct dereferentiation of underlying ptr
    auto const * operator->() const { return &(*link_); }
    friend std::ostream & operator<<(std::ostream & out, LinkPtr const & l) {
        out << *(l.link_);
        return out;
    }
    auto position(size_t i) const { return link_->position(i); }
    auto reverse(size_t i) const { return link_->reverse(i); }
    auto sameDiagonal(LinkPtr const & rhs) const { return link_->sameDiagonal(*(rhs.link_)); }
    auto sameSequences(LinkPtr const & rhs) const { return link_->sameSequences(*rhs); }
    auto sequence(size_t i) const { return link_->sequence(i); }
    auto span() const { return link_->span(); }
    auto const & first() const { return link_->first(); }
    auto const & second() const { return link_->second(); }

private:
    //! Stores the link pointer
    std::shared_ptr<Link> link_;
};

//! Implements hashing function for LinkPtr
struct LinkPtrHash {
    size_t operator()(LinkPtr const & link) const {
        return LinkHash{}(*link);
    }
};

//! Implements hashing function for LinkPtr, ignoring span
struct LinkPtrHashIgnoreSpan {
    size_t operator()(LinkPtr const & link) const {
        return LinkHashIgnoreSpan{}(*link);
    }
};

//! Implements operator() to check if two LinkPtr s are equal, ignoring span
struct LinkPtrEqualIgnoreSpan {
    bool operator()(LinkPtr const & lhs, LinkPtr const & rhs) const {
        return LinkEqualIgnoreSpan{}(*lhs, *rhs);
    }
};

#endif // LINK_H
