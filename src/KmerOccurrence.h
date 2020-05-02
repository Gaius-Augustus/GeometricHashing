#ifndef KMEROCCURRENCE_H
#define KMEROCCURRENCE_H

#include <bitset>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>

#include "json/json.hpp"
#include "CustomHashGeneral.h"
#include "FastaCollection.h"
#include "IdentifierMapping.h"
#include "JsonStream.h"
#include "ReverseComplement.h"

//! Stores the occurrence of a k-mer in 8 byte
class KmerOccurrence {
public:
    //! Constructor
    KmerOccurrence(uint8_t genomeID, uint32_t sequenceID, size_t position, bool reverseStrand, std::string const & kmer)
          // initialize bitvector with ..00 gggg bits
        : data_{genomeID} {
        if (genomeID > 0xf) { throw std::runtime_error("[ERROR] -- KmerOccurrence -- Too many input genomes"); }
        if (sequenceID > 0x3ffff) { throw std::runtime_error("[ERROR] -- KmerOccurrence -- Too many input sequences"); }
        if (position > 0xffffffffff) { throw std::runtime_error("[ERROR] -- KmerOccurrence -- Input sequence too long"); }

        // set remaining bits
        if (reverseStrand) { data_.set(4); }
        if (kmer > reverseComplement(kmer)) { data_.set(5); }

        std::bitset<64> sequenceBits(sequenceID);   //      ..00 00ss ssss ssss ssss ssss
        sequenceBits <<= 6;                         // ..00 ssss ssss ssss ssss ss00 0000
        data_ |= sequenceBits;

        std::bitset<64> positionBits(position);     // same principle as above
        positionBits <<= 24;
        data_ |= positionBits;
    }
    //! KmerOccurrence stores the first position of a k-mer, use this to calculate the central position
    /*! If the k-mer length is even, the center is the left/smaller of both possibilities */
    static size_t centerPosition(size_t firstPosition, size_t k) {
        return firstPosition + static_cast<size_t>(std::ceil(static_cast<double>(k)/2.)) - 1;
    }
    //! getter for data member
    auto const & data() const { return data_; }
    //! Differs from KmerOccurrence::operator==() as k-mer string is ignored
    static bool equalSpot(KmerOccurrence const & lhs, KmerOccurrence const & rhs) {
        return lhs.genome() == rhs.genome()
                && lhs.sequence() == rhs.sequence()
                && lhs.reverse() == rhs.reverse()
                && lhs.position() == rhs.position();
    }
    //! Getter for the genome ID
    uint8_t genome() const { return genomeID_(); }
    //! Get the k-mer string from this occurrence
    std::string getKmerString(FastaCollection const & fastas,
                              IdentifierMapping const & idMap,
                              size_t k) const {
        auto & fastaSequence = fastas.fastaSequence(sequenceID_(), idMap);
        auto kmer = fastaSequence.sequence().substr(this->position(), k);
        return storedKmer(kmer);
    }
    //! Implements operator== by checking if all bits are equal
    bool operator==(KmerOccurrence const & rhs) const {
        return data_ == rhs.data_;
    }
    //! Implements operator< by comparing genomeID, sequenceID, position and strand (in that order)
    bool operator<(KmerOccurrence const & rhs) const {
        // tuples have the desired comparison logic
        std::tuple<uint8_t, uint32_t, size_t, bool> lhsTuple{genome(), sequence(),
                                                              position(), reverse()};
        std::tuple<uint8_t, uint32_t, size_t, bool> rhsTuple{rhs.genome(), rhs.sequence(),
                                                              rhs.position(), rhs.reverse()};
        return lhsTuple < rhsTuple;
    }
    //! Implements operator<< for a KmerOccurrence object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, KmerOccurrence const & occ) {
        out << "(" << static_cast<uint>(occ.genome()) << ", " << occ.sequence() << ", "
            << occ.position() << ", " << occ.reverse() << ")";
        return out;
    }
    //! Getter for the position
    size_t position() const { return position_(); }
    //! Getter for strand information
    bool reverse() const { return reverseStrand_(); }
    //! Getter for the sequence ID
    uint32_t sequence() const { return sequenceID_(); }
    //! Returns the true k-mer at this occurrence
    std::string storedKmer(std::string const & queryKmer) const {
        auto rc = reverseComplement(queryKmer);
        return (data_.test(5))
                ? std::max(queryKmer, rc)
                : std::min(queryKmer, rc);
    }
    //! Return a JsonValue representation
    auto toJsonValue(size_t span, std::shared_ptr<IdentifierMapping const> idMap) const {
        std::array<std::string, 4> occArray{
            std::to_string(centerPosition(position_(), span)),
            std::to_string(reverseStrand_()),
            idMap->queryGenomeName(genomeID_()),
            idMap->querySequenceName(sequenceID_())
        };
        return JsonValue{occArray};
    }

private:
    //! Getter for Genome ID
    uint8_t genomeID_() const {
        std::bitset<64> genomeMask{0xf}; // ..00 1111
        auto idBits = data_ & genomeMask;
        return idBits.to_ulong();
    }
    //! Getter for Position (e.g. absolute position or tile ID)
    size_t position_() const {
        std::bitset<64> genomeMask{0xffffffffff};   // ..00 (1111) (eleven times)
        auto positionBits = data_ >> 24;            // shift s, s, r and g bits away
        positionBits &= genomeMask;
        return positionBits.to_ullong();
    }
    //! Getter for Strand information
    bool reverseStrand_() const {
        return data_.test(4);   // position from right to left
    }
    //! Getter for Sequence ID
    uint32_t sequenceID_() const {
        std::bitset<64> genomeMask{0x3ffff}; // ..00 0011 1111 1111 1111 1111
        auto idBits = data_ >> 6;           // shift 'cr gggg' bits away
        idBits &= genomeMask;
        return idBits.to_ulong();
    }

    //! Stores the data in a bitfield of size 64 to avoid memory waste via padding
    /*! Bitset layout:
     * pppp pppp pppp pppp pppp pppp pppp pppp
     * pppp pppp ssss ssss ssss ssss sscr gggg
     * where \c p are the bits for position, \c s are the bits for sequence id,
     * \c c is a bit that indicates if the lexicographically bigger of kmer and reverse complement
     * is stored (both are treated as the same k-mer and need to be separated this way),
     * \c r is the bit for reverse strand information and \c g are the bits for genome id
     * This implies that this class represent at most
     *   16 different genomes
     *   262,143 different sequences (of all genomes!)
     *   Sequence lengths or tile IDs up to ‭1,099,511,627,775‬ */
    std::bitset<64> data_;
};



//! Computes hash value of a KmerOccurrence object
struct KmerOccurrenceHash {
    size_t operator()(KmerOccurrence const & occ) const {
        std::hash<std::bitset<64>> hashfun;
        return hashfun(occ.data());
    }
};



//! Computes hash value of a KmerOccurrence object ignoring the k-mer bit
struct KmerOccurrencePositionHash {
    size_t operator()(KmerOccurrence const & occ) const {
        std::hash<std::bitset<64>> hashfun;
        auto seed = hashfun(occ.genome());
        customCombineHash(seed, hashfun(occ.sequence()));
        customCombineHash(seed, hashfun(occ.reverse()));
        customCombineHash(seed, hashfun(occ.position()));
        return seed;
    }
};



//! Stores two KmerOccurrence s and provides methods to evaluate their distance
class KmerOccurrencePair {
public:
    //! Default c'tor for parallel sort with tbb
    KmerOccurrencePair() : first_{0,0,0,false,""}, second_{0,0,0,false,""} {}
    //! c'tor
    /*! \param i First KmerOccurrence
     * \param j Second KmerOccurrence
     *
     * \details Stores the lower of \c i and \c j as \c first_ */
    KmerOccurrencePair(KmerOccurrence const & i,
                       KmerOccurrence const & j)
        : first_((i < j) ? i : j), second_((i < j) ? j : i) {}
    //! Returns distance j - i (where i <= j)
    long distance() const {
        auto pos1 = first_.position();
        auto pos2 = second_.position();
        return pos2 - pos1;
    }
    //! Getter for member \c first_
    KmerOccurrence const & first() const { return first_; }
    //! Compares two KmerOccurrenceDistance s (lhs < rhs)
    /*! If \c rhs.first_ is smaller than \c rhs.first, returns true.
     * If both \c first_ s are at the same spot, compare the \c second s */
    bool operator<(KmerOccurrencePair const & rhs) const {
        // Need an interweaved comparison of first and second to ensure that diagonals get sorted together:
        // if this->first.genome < rhs.first.genome
        //   if first genomes same and this->second.genome < rhs.second.genome --> brings all same-genome-matches together
        //
        //     if both genomes same and this->first.sequence < rhs.first.sequence
        //       if both genomes and first sequence same and this->second.sequence < rhs.second.sequence --> brings all same-sequence-matches together
        //
        //         if both genomes and sequences same, compare distance (here we bring diagonals together)
        //           rest by normal comparison

        return ( (first_.genome() < rhs.first_.genome())
                 || ((first_.genome() == rhs.first_.genome()) && (second_.genome() < rhs.second_.genome()))
                 || ((first_.genome() == rhs.first_.genome()) && (second_.genome() == rhs.second_.genome())
                     && (first_.sequence() < rhs.first_.sequence()))
                 || ((first_.genome() == rhs.first_.genome()) && (second_.genome() == rhs.second_.genome())
                     && (first_.sequence() == rhs.first_.sequence()) && (second_.sequence() < rhs.second_.sequence()))
                 || ((first_.genome() == rhs.first_.genome()) && (second_.genome() == rhs.second_.genome())
                     && (first_.sequence() == rhs.first_.sequence()) && (second_.sequence() == rhs.second_.sequence())
                     && (distance() < rhs.distance()))
                 || ((first_.genome() == rhs.first_.genome()) && (second_.genome() == rhs.second_.genome())
                     && (first_.sequence() == rhs.first_.sequence()) && (second_.sequence() == rhs.second_.sequence())
                     && (distance() == rhs.distance())
                     && (first_ < rhs.first_))
                 || ((first_.genome() == rhs.first_.genome()) && (second_.genome() == rhs.second_.genome())
                     && (first_.sequence() == rhs.first_.sequence()) && (second_.sequence() == rhs.second_.sequence())
                     && (distance() == rhs.distance())
                     && (KmerOccurrence::equalSpot(first_, rhs.first_)) && (second_ < rhs.second_)) );
    }
    //! Implements operator<< for a KmerOccurrenceDistance object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, KmerOccurrencePair const & dist) {
        out << dist.first() << " -- " << dist.second() << std::endl;
        return out;
    }
    //! This checks for equality using equalSpot as spaced-seed-strings from the same spot may differ
    bool operator==(KmerOccurrencePair const & rhs) const {
        return KmerOccurrence::equalSpot(first_, rhs.first_)
                && KmerOccurrence::equalSpot(second_, rhs.second_);
    }
    //! \c true if both KmerOccurrence s are on the same diagonal
    bool sameDistance(KmerOccurrencePair const & rhs) const {
        return first_.genome() == rhs.first_.genome()          // same genomes
                && second_.genome() == rhs.second_.genome()
                && first_.sequence() == rhs.first_.sequence()  // same sequences
                && second_.sequence() == rhs.second_.sequence()
                && first_.reverse() == rhs.first_.reverse()    // same strands
                && second_.reverse() == rhs.second_.reverse()
                && distance() == rhs.distance();                // same distance
    }
    //! Getter for member \c second_
    KmerOccurrence const & second() const { return second_; }

protected:
    //! Smaller of both KmerOccurrence s
    KmerOccurrence first_;  // always "smaller" (or equal)
    //! Bigger of both KmerOccurrence s
    KmerOccurrence second_; // always "bigger" (or equal)
};



struct KmerOccurrencePairHash {
    size_t operator()(KmerOccurrencePair const & d) const {
        KmerOccurrencePositionHash hashfun;
        auto seed = hashfun(d.first());
        customCombineHash(seed, hashfun(d.second()));
        return seed;
    }
};

#endif // KMEROCCURRENCE_H
