#ifndef TWOBITKMER_H
#define TWOBITKMER_H

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <functional> // std::hash
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <tsl/hopscotch_set.h>
#include "CustomHashGeneral.h"
#include "ReverseComplement.h"



//! Custom exception to specifically catch if a non-ACGTU nucleotide occured in input k-mer to TwoBitKmerHash
class TwoBitKmerWrongBaseException : public std::runtime_error {
public:
    explicit TwoBitKmerWrongBaseException(std::string const & what_arg) : std::runtime_error(what_arg) {}
};


template <typename TwoBitKmerData>
struct TwoBitKmerHash;  // fwd declare

//! Succinct representation of a DNA k-mer with two bits per base
/*! Only valid bases are A, C, G and T (or U) in upper case,
 * every other character leads to an exception.
 * If input contained 'U's, they are lost and converted to 'T' when
 * recovered from this representation.
 *
 * On typical modern 64 bit CPUs, due to padding the lowest allocation unit is
 * often(*) 8 Byte, even if space for a single Byte (e.g. std::bitset<8>) is queried.
 * Thus, this class works with multiples of std::bitset<64>, i.e. (up to) 32 characters (bases).
 *
 * TwoBitKmerData classes need to implement the methods
 * * TwoBitKmerData(size_t size)
 * * std::bitset<64> & bitset(size_t position) // return writable bitset ref corresponding to \c position
 * * std::bitset<64> const & bitsetRO(size_t position) // read-only reference
 * * bool operator==(TwoBitKmerData const & rhs) const
 * * size_t objectSize() const
 * * size_t size() const
 *
 * (*)explanation: http://www.catb.org/esr/structure-packing/ */
template<typename TwoBitKmerData>
class TwoBitKmer {
public:
    TwoBitKmer(std::string const & kmer)
        : bitset_{kmer.size()} {
        for (size_t position = 0; position < kmer.size(); ++position) {
            setBase(bitset_.bitset(position), position, kmer[position]);
        }
    }
    TwoBitKmer(std::string const & kmer, std::shared_ptr<bool> const & validFlag)
        : bitset_{kmer.size()} {
        for (size_t position = 0; position < kmer.size(); ++position) {
            if ((kmer[position] != 'A') && (kmer[position] != 'C') && (kmer[position] != 'G') && (kmer[position] != 'T') && (kmer[position] != 'U')) {
                *validFlag = false;
                return;
            }
            setBase(bitset_.bitset(position), position, kmer[position]);
        }
        *validFlag = true;
    }
    //! Alias for base()
    char at(size_t position) const { return base(position); }
    //! Return the base character at \c position
    char base(size_t position) const {
        return getBase(bitset_.bitsetRO(position), position);
    }
    //! Return the length k of the k-mer
    size_t length() const {
        return bitset_.size();
    }
    //! Calculate memory consumption of this object in bytes
    size_t objectSize() const {
        return bitset_.objectSize();
    }
    //! Compare for equality
    bool operator==(TwoBitKmer<TwoBitKmerData> const & rhs) const { return bitset_ == rhs.bitset_; }
    //! Compare for equality with a \c std::string
    bool operator==(std::string const & rhs) const { return toString() == rhs; }
    //! Check if this is less than another TwoBitKmer instance
    bool operator<(TwoBitKmer const & rhs) const {
        return toString() < rhs.toString();
    }
    //! Create TwoBitKmer of \c this reverse complement
    TwoBitKmer reverseComplement() const {
        return TwoBitKmer( ::reverseComplement(toString()) );
    }
    //! Create k-mer string from 2bit representation
    std::string toString() const {
        std::string kmer;
        for (size_t i = 0; i < length(); ++i) {
            kmer += base(i);
        }
        return kmer;
    }
    //! Convenient \c iostream output
    friend std::ostream & operator<<(std::ostream & out, TwoBitKmer const & tbk) {
        out << tbk.toString();
        return out;
    }

protected:
    friend struct TwoBitKmerHash<TwoBitKmerData>;

    //! Amount of bits to shift a single \c std::bitset<64> to get 7*(0000) 00xx where xx are the resp. bits for the base at \c position
    size_t bitsetShift(size_t position) const {
        /* base sequence in bitvector from right to left, i.e. '... 3322 1100'
         * thus, shift twoBitCode accordingly:
         * first base: no shift, 0000 00xx
         * second base: << 2 --> 0000 xx00
         * third base:  << 4 --> 00xx 0000
         * and so on */
        return 2 * (position % 32);
    }
    //! Return the base character at \c position
    char getBase(std::bitset<64> bitset, size_t position) const {
        auto shift = bitsetShift(position);
        bitset >>= shift;   // move relevant bits to rightmost bit-pair
        std::bitset<64> outputMask(3);  // 7*(0000) 0011
        bitset &= outputMask;           // 7*(0000) 00xx
        switch(bitset.to_ulong()) {
            case 0 : return 'A';
            case 1 : return 'C';
            case 2 : return 'T';
            case 3 : return 'G';
            default : throw std::runtime_error("[ERROR] -- TwoBitKmer::getBase() -- Something went wrong");
        }
    }
    //! Set the bits to the resp. base enconding of the base at \c position
    void setBase(std::bitset<64> & bitset, size_t position, char base) {
        if ((base != 'A') && (base != 'C') && (base != 'G') && (base != 'T') && (base != 'U')) {
            //throw std::runtime_error("[ERROR] -- TwoBitKmer::setBase() -- Invalid input character");
            throw TwoBitKmerWrongBaseException("[ERROR] -- TwoBitKmer::setBase() -- Invalid input character");
        }
        std::bitset<64> insertMask(3);   // 7*(0000) 0011
        std::bitset<64> baseBits(base);
        /* encoding from ASCII codes: https://stackoverflow.com/a/39244735/5586698
         * A - 0100 0|00|1 --> 7*(0000) 0000
         * C - 0100 0|01|1 --> 7*(0000) 0001
         * G - 0100 0|11|1 --> 7*(0000) 0011
         * T - 0101 0|10|0 --> 7*(0000) 0010
         * U - 0101 0|10|1 --> 7*(0000) 0010 */
        baseBits >>= 1;
        baseBits &= insertMask;
        auto bitShift = bitsetShift(position);
        insertMask <<= bitShift;    // 11 at insert position, 0 otherwise
        baseBits <<= bitShift;      // actual value at insert position, 0 otherwise
        baseBits |= ~insertMask;    // actual value at insert position, 1 otherwise (operator~ flips bits)
        bitset |= insertMask;   // set bits at insert position to 11    (| with 0 leaves other bits unchanged)
        bitset &= baseBits;     // set bits at insert position to value (& with 1 leaves other bits unchanged)
    }

    //! Data member
    TwoBitKmerData bitset_;
};



//! Can hold up to 29-mers in 8 bytes
class TwoBitKmerDataShort {
public:
    TwoBitKmerDataShort(size_t size)
        : bitset_{} {
        if (size > 29) { throw std::runtime_error("[ERROR] -- TwoBitKmerDataShort -- You must use TwoBitKmerDataMedium or TwoBitKmerDataLong for k > 29"); }
        auto sizeBits = std::bitset<64>{size};  // 0000 0000 .... 00xx xxxx
        sizeBits <<= 58;                        // xxxx xx00 .... 0000 0000
        bitset_ |= sizeBits;
    }
    std::bitset<64> & bitset(size_t position) {
        if (position > 29) { throw std::runtime_error("[ERROR] -- TwoBitKmerDataShort::bitset -- You must use TwoBitKmerDataMedium or TwoBitKmerDataLong for k > 29"); }
        return bitset_;
    }
    std::bitset<64> const & bitsetRO(size_t position) const {
        if (position > 29) { throw std::runtime_error("[ERROR] -- TwoBitKmerDataShort::bitsetRO -- You must use TwoBitKmerDataMedium or TwoBitKmerDataLong for k > 29"); }
        return bitset_;
    }
    size_t objectSize() const {
        return sizeof(bitset_);
    }
    bool operator==(TwoBitKmerDataShort const & rhs) const { return bitset_ == rhs.bitset_; }
    size_t size() const {
        auto sizeBits = bitset_;
        sizeBits >>= 58;
        return sizeBits.to_ullong();
    }
private:
    std::bitset<64> bitset_;
};



//! Can hold up to 61-mers in 16 bytes (avoid std::vector overhead)
class TwoBitKmerDataMedium {
public:
    TwoBitKmerDataMedium(size_t size)
        : bitset1_{}, bitset2_{} {
        if (size > 61) { throw std::runtime_error("[ERROR] -- TwoBitKmerDataMedium -- You must use TwoBitKmerDataLong for k > 61"); }
        if (size <= 29) { std::cerr << "[WARNING] -- TwoBitKmerDataMedium -- Use TwoBitKmerDataShort for k <= 29" << std::endl; }
        auto sizeBits = std::bitset<64>{size};  // 0000 0000 .... 00xx xxxx
        sizeBits <<= 58;                        // xxxx xx00 .... 0000 0000
        bitset2_ |= sizeBits;
    }
    std::bitset<64> & bitset(size_t position) {
        return ((position/32) == 0)
                ? bitset1_
                : bitset2_;
    }
    std::bitset<64> const & bitsetRO(size_t position) const {
        return ((position/32) == 0)
                ? bitset1_
                : bitset2_;
    }
    size_t objectSize() const {
        return 2 * sizeof(bitset1_);
    }
    bool operator==(TwoBitKmerDataMedium const & rhs) const {
        return bitset1_ == rhs.bitset1_
                && bitset2_ == rhs.bitset2_;
    }
    size_t size() const {
        auto sizeBits = bitset2_;
        sizeBits >>= 58;
        return sizeBits.to_ullong();
    }
private:
    std::bitset<64> bitset1_;   // bases 1 - 32
    std::bitset<64> bitset2_;   // bases 33 - 61 + size info
};



//! Can hold an arbitrary long k-mer, use only for k > 61
class TwoBitKmerDataLong {
public:
    TwoBitKmerDataLong(size_t size)
        : size_{size},
          // init bitset_ vector with enough bitsets to fit kmer
          bitsetVector_(static_cast<size_t>(std::ceil(static_cast<double>(size)/32.))) {
         if (size <= 61) { std::cerr << "[WARNING] -- TwoBitKmerLong -- You should use TwoBitKmerMedium or TwoBitKmerShort for k <= 61" << std::endl; }
    }
    std::bitset<64> & bitset(size_t position) {
        return bitsetVector_.at(bitsetIndex(position));
    }
    std::bitset<64> const & bitsetRO(size_t position) const {
        return bitsetVector_.at(bitsetIndex(position));
    }
    size_t objectSize() const {
        return sizeof(size_) + (bitsetVector_.size() * sizeof(std::bitset<64>));
    }
    bool operator==(TwoBitKmerDataLong const & rhs) const { return size_ == rhs.size_
                                                                       && bitsetVector_ == rhs.bitsetVector_; }
    size_t size() const {
        return size_;
    }
private:
    //! Index of bitset that holds the resp. bits for the base at \c position
    size_t bitsetIndex(size_t position) const { return position/32; }

    size_t size_;
    std::vector<std::bitset<64>> bitsetVector_;
};



//! Computes hash value of a TwoBitKmer object
template <typename TwoBitKmerData>
struct TwoBitKmerHash {
    size_t operator()(TwoBitKmer<TwoBitKmerData> const & tbk) const {
        std::hash<std::string> stringhash;
        auto hashValue = stringhash(tbk.toString());
        return hashValue;
    }
};



//! Computes a hash that is identical for a TwoBitKmer and its reverse complement
template<typename TwoBitKmerData>
struct TwoBitKmerRCIncludingHash {
    TwoBitKmerHash<TwoBitKmerData> twoBitKmerHash;
    size_t operator()(TwoBitKmer<TwoBitKmerData> const & tbk) const {
        auto tbkRC = tbk.reverseComplement();
        size_t seed = 0;
        if (tbk < tbkRC) {
            customCombineHash(seed, twoBitKmerHash(tbk));
            customCombineHash(seed, twoBitKmerHash(tbkRC));
        } else {
            customCombineHash(seed, twoBitKmerHash(tbkRC));
            customCombineHash(seed, twoBitKmerHash(tbk));
        }
        return seed;
    }
};



//! Compares two TwoBitKmer objects and including both reverse complements
/*! If either of lhs kmer or its respective reverse complement is equal to
 * the rhs kmer or its respective rc, returns true */
template<typename TwoBitKmerData>
struct TwoBitKmerRCIncludingEqual {
    bool operator()(TwoBitKmer<TwoBitKmerData> const & lhs,
                    TwoBitKmer<TwoBitKmerData> const & rhs) const {
        auto lhsRC = lhs.reverseComplement();
        auto rhsRC = rhs.reverseComplement();
        return lhs == rhs
                || lhsRC == rhs
                || lhs == rhsRC     // I think the last two checks are unneccessary...
                || lhsRC == rhsRC;
    }
};



//! Compares two TwoBitKmer objects w.r.t literal order and including both reverse complements
/*! If lhs kmer or its respective reverse complement is lower than the rhs kmer or its
 * respective rc, returns true */
template<typename TwoBitKmerData>
struct TwoBitKmerRCIncludingLess {
    bool operator()(TwoBitKmer<TwoBitKmerData> const & lhs,
                    TwoBitKmer<TwoBitKmerData> const & rhs) const {
        auto lhsRC = lhs.reverseComplement();
        auto rhsRC = rhs.reverseComplement();
        return std::min(lhs, lhsRC) < std::min(rhs, rhsRC);
    }
};

#endif // TWOBITKMER_H
