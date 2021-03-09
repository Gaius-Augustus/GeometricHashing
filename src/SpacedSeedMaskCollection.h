#ifndef SPACEDSEEDMASKCOLLECTION_H
#define SPACEDSEEDMASKCOLLECTION_H

#include <cstdlib>
#include <iostream>
#include <set>
#include <random>
#include <vector>

#include "boost/math/special_functions/binomial.hpp"
#include "prettyprint.hpp"
#include "optimalSpacedSeeds.h"
#include "SpacedSeedMask.h"
#include "StrongType.h"

//! Create and store \c m distinct SpaceSeedMask s
class SpacedSeedMaskCollection {
public:
    // https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
    using Weight = NamedType<size_t, struct WeightTag>;
    using Span = NamedType<size_t, struct SpanTag>;
    using SeedSetSize = NamedType<size_t, struct SeedSetSizeTag>;

    //! c'tor (1)
    /*! \param weight number of set bits
     * \param span size of the bitset
     * \param seedSetSize Number of SpacedSeedMask s to create
     *
     * \details Initializes all members and creates random distinct SpacedSeedMask s */
    SpacedSeedMaskCollection(Weight weight, Span span, SeedSetSize seedSetSize)
        : seedSetSize_{seedSetSize.get()}, spacedSeedMasks_{}, span_{span.get()}, weight_{weight.get()} {
        // max value of m is the binomial coefficient of span over weight
        if (static_cast<double>(seedSetSize_) > boost::math::binomial_coefficient<double>(span_, weight_)) {
            throw std::runtime_error("[ERROR] -- SpacedSeedMaskCollection -- m too large");
        }
        // create m distinct spacedSeedMasks
        std::random_device rd;  // only acquire this once as this may be an expensive system call
        std::default_random_engine rng(rd());
        std::set<SpacedSeedMask> seen;
        for (size_t i = 0; i < seedSetSize_; ++i) {
            bool distinct = false;
            while(!distinct) {
                auto mask = SpacedSeedMask(weight_, span_, rng);
                if (seen.find(mask) == seen.end()) {
                    distinct = true;
                    seen.insert(mask);
                    spacedSeedMasks_.push_back(mask);
                }
            }
        }
    }
    //! c'tor (2)
    /*! \param strv vector of string patterns to create bitsets from
     *
     * \details Initializes all members and creates SpacedSeedMask s according to \c strv */
    SpacedSeedMaskCollection(std::vector<std::string> const & strv)
        : seedSetSize_{strv.size()}, spacedSeedMasks_{}, span_{}, weight_{} {
        for (auto&& str : strv) {
            spacedSeedMasks_.emplace_back(str);
        }
        weight_ = (seedSetSize_ > 0) ? spacedSeedMasks_.at(0).weight() : 0;
        span_ = 0;
        for (auto&& mask : spacedSeedMasks_) {
            if (mask.weight() != weight_) {
                throw std::runtime_error("[ERROR] -- SpacedSeedMaskCollection -- patterns must all have the same weight");
            }
            span_ = (mask.span() > span_) ? mask.span() : span_;    // store biggest span
        }
    }
    //! c'tor (3)
    /*! \param weight Weight of optimal seed
     *
     * \details Tries to load pre-computed optimal seed (span depends on seed).
     * Throws if no pre-computed seed exists for that weight */
    SpacedSeedMaskCollection(Weight weight, SeedSetSize size = SeedSetSize(1))
        : seedSetSize_{optimalSpacedSeeds(weight.get(), size.get()).size()}, spacedSeedMasks_{optimalSpacedSeeds(weight.get(), size.get())},
          span_{0}, weight_{weight.get()} {
        for (auto&& mask : spacedSeedMasks_) {
            span_ = (mask.span() > span_) ? mask.span() : span_;    // store biggest span
        }
    }

    //! Getter for member \c spacedSeedMasks_
    auto const & masks() const { return spacedSeedMasks_; }
    //! Transform \c spacedSeedMasks_ to a vector of string representations
    auto masksAsString() const {
        std::vector<std::string> maskStr;
        for (auto&& mask : spacedSeedMasks_) {
            maskStr.emplace_back(mask.string());
        }
        return maskStr;
    }
    //! Getter for member \c span_
    auto maxSpan() const { return span_; }
    //! Getter for member \c m_
    auto size() const { return seedSetSize_; }
    //! Forward getter for span of seed \c i
    auto span(size_t i) const { return spacedSeedMasks_.at(i).span(); }
    //! Getter for member \c l_
    auto weight() const { return weight_; }

    friend std::ostream & operator<<(std::ostream & out, SpacedSeedMaskCollection const & c) {
        out << c.spacedSeedMasks_;
        return out;
    }

private:
    //! Number of SpacedSeedMask s
    size_t seedSetSize_;
    //! Vector of all SpacedSeedMask s
    std::vector<SpacedSeedMask> spacedSeedMasks_;
    //! Biggest span of the (spaced) seed(s)
    size_t span_;
    //! Number of match positions (1s) in the seed(s)
    size_t weight_;
};

#endif // SPACEDSEEDMASKCOLLECTION_H
