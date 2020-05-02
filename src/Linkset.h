#ifndef LINKSET_H
#define LINKSET_H

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <filesystem>
#include <fstream>  // for writing statistics to file
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <unordered_set>
#include <vector>

#include "hopscotch-map/hopscotch_map.h"
#include "json/json.hpp"
#include "prettyprint/prettyprint.hpp"
#include "IdentifierMapping.h"
#include "Link.h"
#include "MemoryMonitor.h"
#include "ProgressBar.h"
#include "Timestep.h"

//! Class that interacts with the underlying graph and extracts Link s from it
/*! Loads the graph into memory (in fact, calls the methods provided by the
 * interface to get the graph loaded into memory), calls the methods from
 * the interface to iterate over all k-mers (with the right properties).
 * If a k-mer occurs in the reference genome and in at least one other
 * genome, all Link s between the k-mer occurences are created an stored
 * together with the Link counts */
class Linkset {
public:
    /* using LinksetType = std::unordered_map<std::shared_ptr<Link const>, size_t,
                                           LinkPtrHash, LinkPtrEqual>; */
    using LinksetType = tsl::hopscotch_map<std::shared_ptr<Link const>, size_t,
                                           LinkPtrHash, LinkPtrEqual>;
    //! Constructor (1)
    /*! \param indetifierMapping Instance of IdentifierMapping that already
     * knows about all IDs
     * \param linkLimit If a k-mer emits more than \c linkLimit Link s, store only a
     * random Link subset of size \c linkLimit
     * \param occurrencePerGenomeMax If a k-mer occurs more than \c occurrencePerGenomeMax times in any genome,
     * discard it completely
     * \param occurencePerGenomeMin If a k-mer occurs less than \c occurencePerGenomeMin times in any genome,
     * discard it completely
     * \param keepKmers If \c true , create HeavyLink s that store the respective k-mers instead of lightweight Link s
     *
     * \details Creates an empty Linkset with a user-defined IdentifierMapping
     * member that must already know all names of input genomes and sequences. */
    explicit Linkset(std::shared_ptr<IdentifierMapping const> identifierMapping,
                     size_t linkLimit, bool linkLimitDiscard, size_t occurrencePerGenomeMax,
                     size_t occurrencePerGenomeMin)
        : discardedNotInReference_{0}, discardedOnlyInReference_{0},
          discardedTooFewOccurrences_{0}, discardedTooManyOccurrences_{0},
          discardExceeding_{linkLimitDiscard},
          idMapping_{identifierMapping}, linkLimit_{linkLimit},
          linkset_{}, occurrencePerGenomeMax_{occurrencePerGenomeMax},
          occurrencePerGenomeMin_{occurrencePerGenomeMin}, rd_{} {}
    //! Add a Link into the Linkset or increase counter for this Link
    void addLink(std::shared_ptr<Link const> link);
    //! Create a Link in the Linkset from a vector of occurrences
    void createLinks(std::vector<KmerOccurrence> const & occurrences);
    //! Getter for the member variable \c idMapping_
    std::shared_ptr<IdentifierMapping const> idMapping() const { return idMapping_; }
    //! Return number of times that \c link was created
    size_t linkCount(std::shared_ptr<Link const> const & link) const {
        return (linkset_.find(link) != linkset_.end())
                ? linkset_.at(link)
                : 0;
    }
    //! Getter for the member variable \c linkset_
    auto const & linkset() const { return linkset_; }
    //! Return number of discarded k-mers during construction of this Linkset
    size_t numDiscardedKmers() const { return discardedNotInReference_ + discardedOnlyInReference_ + discardedTooFewOccurrences_ + discardedTooManyOccurrences_; }
    //! Return number of genomes in input data
    size_t numGenomes() const { return idMapping_->numGenomes(); }
    //! Return number of Link s in this Linkset
    size_t numLinks() const { return linkset_.size(); }
    Linkset& operator=(Linkset && other) {
        idMapping_ = other.idMapping_;
        linkLimit_ = std::move(other.linkLimit_);
        occurrencePerGenomeMax_ = std::move(other.occurrencePerGenomeMax_);
        occurrencePerGenomeMin_ = std::move(other.occurrencePerGenomeMin_);
        linkset_ = std::move(other.linkset_);
        return *this;
    }
    //! Implements operator== for Linkset
    /*! Checks if \c linkset_ members are equal, i.e. the same Link s with the same counts
     * must be present, plus the \c idMapping_ members and the remaining members must be equal as well */
    bool operator==(Linkset const & rhs) const {    // to test if load/save works correctly
        for (auto&& elem : linkset_) {  // unordered_map::operator==() does not work correctly because of shared_ptrs
            if (rhs.linkset_.find(elem.first) == rhs.linkset_.end()
                    || elem.second != rhs.linkset_.at(elem.first)) {
                return false;
            }
        }
        for (auto&& elem : rhs.linkset_) {
            if (linkset_.find(elem.first) == linkset_.end()
                    || elem.second != linkset_.at(elem.first)) {
                return false;
            }
        }
        return *idMapping_ == *(rhs.idMapping_)
                && linkLimit_ == rhs.linkLimit_
                && occurrencePerGenomeMax_ == rhs.occurrencePerGenomeMax_
                && occurrencePerGenomeMin_ == rhs.occurrencePerGenomeMin_;
    }
    //! Implements operator<< for Linkset objects for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, Linkset const & ls) {
        out << "ID Mapping:" << std::endl << *(ls.idMapping_) << std::endl;
        out << "Link set:" << std::endl;
        for (auto&& element : ls.linkset_) {
            out << "Count" << element.second << std::endl;
            out << *(element.first) << std::endl;
        }
        return out;
    }

private:
    //! Helper method that constructs all Link s from a k-mer
    void addLinksFromKmer(std::vector<KmerOccurrence> const & occurrences);

    //! Count discarded or skipped k-mers (not on reference)
    size_t discardedNotInReference_;
    //! Count discarded or skipped k-mers (only on reference)
    size_t discardedOnlyInReference_;
    //! Count discarded or skipped k-mers (not only on reference but not right number of occurences)
    size_t discardedTooFewOccurrences_;
    //! Count discarded or skipped k-mers (not only on reference but not right number of occurences)
    size_t discardedTooManyOccurrences_;
    //! Discard seeds with too many Links rather than sampling
    bool discardExceeding_;
    //! Assigns IDs to genome and sequence strings in order of their appeareances
    /*! The reference genome always gets ID '0' */
    std::shared_ptr<IdentifierMapping const> idMapping_;
    //! Maximum number of links per k-mer to be created
    size_t linkLimit_;
    //! Stores a mapping from Link s to their counts
    LinksetType linkset_;
    //! Maximum number of k-mer occurences per genome
    size_t occurrencePerGenomeMax_;
    //! Minimum number of k-mer occurences per genome
    size_t occurrencePerGenomeMin_;
    //! Used to obtain seed for random link selection if there are too many possibilities
    std::random_device rd_;
};

#endif // LINKSET_H
