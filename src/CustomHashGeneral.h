#ifndef CUSTOMHASHGENERAL_H
#define CUSTOMHASHGENERAL_H

#include <array>
#include <cstdlib>
#include <vector>



//! Copied from boost 1.69.0 function 'void hash_combine_impl(SizeT& seed, SizeT value)' from boost/functional/hash/hash.hpp
/*! Boost does not allow direct passing of hash values in hash_combine from the user interface,
 * so just DIY here */
inline void customCombineHash(size_t & seed, size_t hash) {
    seed ^= hash + 0x9e3779b9 + (seed<<6) + (seed>>2);
}



//! Computes Hash value of an array
template <typename T, size_t N, typename THash>
struct ArrayHash {
    size_t operator()(std::array<T const, N> const & a) const {
        size_t seed = 0;
        THash hashfun;
        for (size_t i = 0; i < N; ++i) {
            customCombineHash(seed, hashfun(a[i]));
        }
        return seed;
    }
};


//! Computes Hash value of a vector
struct VectorHash {
    template <typename T>
    size_t operator()(std::vector<T> const & v) const {
        std::hash<T> hashfun;
        size_t seed = 0;
        for (auto&& elem : v) {
            customCombineHash(seed, hashfun(elem));
        }
        return seed;
    }
};

#endif // CUSTOMHASHGENERAL_H
