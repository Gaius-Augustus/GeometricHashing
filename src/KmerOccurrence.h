#ifndef KMEROCCURRENCE_H
#define KMEROCCURRENCE_H

#include <bitset>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>

#include "mabl3/JsonStream.h"
#include "nlohmann/json.hpp"
#include "CustomHashGeneral.h"
#include "FastaCollection.h"
#include "IdentifierMapping.h"
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
        std::tuple<uint8_t, uint32_t, size_t, bool, bool> lhsTuple{genome(), sequence(),
                                                                   position(), reverse(),
                                                                   data_.test(5)};
        std::tuple<uint8_t, uint32_t, size_t, bool, bool> rhsTuple{rhs.genome(), rhs.sequence(),
                                                                   rhs.position(), rhs.reverse(),
                                                                   rhs.data_.test(5)};
        return lhsTuple < rhsTuple;
    }
    //! Implements operator<< for a KmerOccurrence object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, KmerOccurrence const & occ) {
        out << "(" << static_cast<uint>(occ.genome()) << ", " << occ.sequence() << ", "
            << occ.position() << ", " << occ.reverse() << "){" << occ.data_.test(5) << "}";
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
    auto toJsonValue(size_t span, IdentifierMapping const & idMap) const {
        std::array<std::string, 5> occArray{
            std::to_string(position_()),//std::to_string(centerPosition(position_(), span)),
            std::to_string(reverseStrand_()),
            idMap.queryGenomeName(genomeID_()),
            idMap.querySequenceName(sequenceID_()),
            std::to_string(span)
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



//! To tell containers to use equalSpot instead of operator== with KmerOccurrence s
struct KmerOccurrencePositionEqual {
    bool operator()(KmerOccurrence const & lhs, KmerOccurrence const & rhs) const {
        return KmerOccurrence::equalSpot(lhs, rhs);
    }
};

#endif // KMEROCCURRENCE_H
