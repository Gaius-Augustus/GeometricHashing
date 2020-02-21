#ifndef SPACEDSEEDMASK_H
#define SPACEDSEEDMASK_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#include "boost/dynamic_bitset.hpp"

//! Basically a wrapper around \c boost::dynamic_bitset that creates such a bitset with \c weight randomly set bits
class SpacedSeedMask {
public:    
    //! C'tor (1)
    /*! \param weight number of set bits
     * \param span size of the bitset
     *
     * \details Acquires a std::random_device, creates a \c boost::dynamic_bitset and randomly sets \c weight of the bits to 1.
     * Throws if \c weight is larger than \c span */
    SpacedSeedMask(size_t weight, size_t span)
        : mask_{span} {
        if (weight > span) { throw std::runtime_error("[ERROR] -- SpacedSeedMask -- weight > span"); }
        if (weight == 0) { throw std::runtime_error("[ERROR] -- SpacedSeedMask -- weight == 0"); }
        std::random_device rd;
        std::default_random_engine rng(rd());
        createRandomMask(weight, span, rng);
    }
    //! C'tor (3)
    /*! \param weight number of set bits
     * \param span size of the bitset
     * \param rng reference to a rng to use instead of acquiring a new one
     *
     * \details Creates a \c boost::dynamic_bitset and randomly sets \c weight of the bits to 1.
     * Throws if \c weight is larger than \c span */
    SpacedSeedMask(size_t weight, size_t span, std::default_random_engine & rng)
        : mask_{span} {
        if (weight > span) { throw std::runtime_error("[ERROR] -- SpacedSeedMask -- weight > span"); }
        if (weight == 0) { throw std::runtime_error("[ERROR] -- SpacedSeedMask -- weight == 0"); }
        createRandomMask(weight, span, rng);
    }
    //! C'tor (3)
    /*! \param str string of 0s and 1s to construct bitset from
     *
     * \details Creates a \c boost::dynamic_bitset and sets the bit according to \c str
     * If str is '1001001100', setPositions() returns {0, 3, 6, 7}, i.e. from left-to-right */
    SpacedSeedMask(std::string const & str)
        : mask_{std::string(str.rbegin(), str.rend())} {}
    //! Return span by querying from \c mask_
    auto span() const { return mask_.size(); }
    //! Return string representation
    auto string() const {
        std::string s{""};
        for (size_t i = 0; i < mask_.size(); ++i) {
            s += (mask_.test(i)) ? "1" : "0";
        }
        return s;
    }
    //! Return weight by querying from \c mask_
    auto weight() const { return mask_.count(); }
    //! Return const reference to \c mask_
    auto const & mask() const { return mask_; }
    //! Equality comparison
    bool operator==(SpacedSeedMask const & rhs) const { return mask_ == rhs.mask_; }
    //! Less comparison
    bool operator<(SpacedSeedMask const & rhs) const { return mask_ < rhs.mask_; }
    //! Returns a vector of positions of the set bits
    auto getSetPositions() const {
        std::vector<size_t> pos;
        for (size_t i = 0; i < mask_.size(); ++i) {
            if (mask_.test(i)) { pos.emplace_back(i); }
        }
        return std::move(pos);
    }
    //! Fwd \c test() method of \c mask_ for bit at \c pos
    auto test(size_t pos) const { return mask_.test(pos); }

    friend std::ostream & operator<<(std::ostream & out, SpacedSeedMask const & m) {
        out << m.mask_;
        return out;
    }

private:
    //! Factory function to create the mask
    inline void createRandomMask(size_t weight, size_t span, std::default_random_engine & rng) {
        std::uniform_int_distribution<size_t> unif(0, (span-1));
        while (mask_.count() < weight) {
            auto n = unif(rng);
            mask_.set(n);
        }
    }

    //! Bitmask representing the spaced seed
    boost::dynamic_bitset<> mask_;
};

#endif // SPACEDSEEDMASK_H
