#ifndef IDENTIFIERMAPPING_H
#define IDENTIFIERMAPPING_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "prettyprint.hpp"
#include <tsl/hopscotch_map.h>
#include "CustomHashGeneral.h"



//! Stores a mapping of genome and sequence IDs to the identifier strings
/*! The reference genome always gets genome ID 0. The other IDs are assigned
 * by incrementing the last ID each time a new string is queried */
class IdentifierMapping {
public:
    struct SequenceGenomeTuple {
        size_t gid;
        std::string sequence;
        bool operator==(SequenceGenomeTuple const & rhs) const {
            return gid == rhs.gid && sequence == rhs.sequence;
        }
        bool operator!=(SequenceGenomeTuple const & rhs) const {
            return !(*this == rhs);
        }
        friend std::ostream & operator<<(std::ostream & out, SequenceGenomeTuple const & t) {
            out << "(" << t.gid << ", " << t.sequence << ")";
            return out;
        }
    };
    struct SequenceGenomeTupleHash {
        size_t operator()(SequenceGenomeTuple const & t) const {
            size_t seed = 0;
            customCombineHash(seed, std::hash<size_t>{}(t.gid));
            customCombineHash(seed, std::hash<std::string>{}(t.sequence));
            return seed;
        }
    };

    //! Constructor (1)
    /*! \details creates an empty idMap */
    IdentifierMapping()
        : genomeIDToName_{},
          genomeNameToID_{},
          nextGenomeID_{0},
          nextSequenceID_{0},
          sequenceIDToTuple_{},
          sequenceTupleToID_{} {}
    //! Constructor (2)
    /*! \param referenceGenome Identifier string of the reference genome
     *
     * \details Assigns genome ID 0 to the reference genome */
    IdentifierMapping(std::string const & referenceGenome)
        : genomeIDToName_{},
          genomeNameToID_{},
          nextGenomeID_{0},
          nextSequenceID_{0},
          sequenceIDToTuple_{},
          sequenceTupleToID_{}  {
        // reference genome always gets ID '0'
        insertNewGenome(referenceGenome);
    }
    //! Copies state of the source IdentifierMapping
    void completeUpdate(IdentifierMapping const & source) {
        *this = source;
    }
    //! Const getter for genomeIDToName_
    auto const & genomeIDToName() const { return genomeIDToName_; }
    //! return true if genome has an ID assigned
    bool genomeKnown(std::string const & genome) const {
        return genomeNameToID_.find(genome) != genomeNameToID_.end();
    }
    //! Const getter for genomeNameToID_
    auto const & genomeNameToID() const { return genomeNameToID_; }
    //! Returns the total number of genomes in the input data
    size_t numGenomes() const { return nextGenomeID_; }
    //! Returns the total number of sequences in the input data
    size_t numSequences() const { return nextSequenceID_; }
    //! Memberwise comparison
    bool operator==(IdentifierMapping const & rhs) const {
        return genomeIDToName_ == rhs.genomeIDToName_
                && genomeNameToID_ == rhs.genomeNameToID_
                && nextGenomeID_ == rhs.nextGenomeID_
                && nextSequenceID_ == rhs.nextSequenceID_
                && sequenceIDToTuple_ == rhs.sequenceIDToTuple_
                && sequenceTupleToID_ == rhs.sequenceTupleToID_;
    }
    //! Returns the ID of a genome (1)
    /*! If the query is not present in the map, it is inserted now
     * Increments reference query count if parameter is query identifier string */
    size_t queryGenomeID(std::string const & genome) {
        insertNewGenome(genome);
        return genomeNameToID_.at(genome);
    }
    //! Returns the ID of a genome (2) -- throws it genome not known
    size_t queryGenomeIDConst(std::string const & genome) const {
        try {
            return genomeNameToID_.at(genome);
        } catch (std::exception const & e) {
            std::cerr << "[ERROR] -- IdentifierMapping::queryGenomeIDConst -- genome " + genome + " not found" << std::endl;
            throw e;
        }
    }
    size_t queryGenomeIDFromSequence(size_t sid) const {
        try {
            return sequenceIDToTuple_.at(sid).gid;
        } catch (std::exception const & e) {
            std::cerr << "[ERROR] -- IdentifierMapping::queryGenomeIDFromSequence -- sequence " + std::to_string(sid) + " not found" << std::endl;
            throw e;
        }
    }
    //! Returns the identifier string of a genome
    /*! If the query is not present in the vector, it throws */
    std::string const & queryGenomeName(size_t const id) const {
        try {
            return genomeIDToName_.at(id);   // throws if ID not present
        } catch (std::exception const & e) {
            std::cerr << "[ERROR] -- IdentifierMapping::queryGenomeName -- genome " + std::to_string(id) + " not found" << std::endl;
            throw e;
        }
    }
    //! Returns the ID of a sequence (1)
    /*! If the query is not present in the map, it is inserted now */
    size_t querySequenceID(std::string const & sequence, std::string const & genome) {
        insertNewGenome(genome);    // make sure genome has a gid
        auto gid = genomeNameToID_.at(genome);
        insertNewSequence(sequence, gid);   // inserts sequence if unknown
        return sequenceTupleToID_.at(SequenceGenomeTuple{gid, sequence});
    }
    //! Returns the ID of a sequence (2)
    size_t querySequenceIDConst(std::string const & sequence, std::string const & genome) const {
        try {
            auto gid = genomeNameToID_.at(genome);
            return sequenceTupleToID_.at(SequenceGenomeTuple{gid, sequence});
        } catch (std::exception const & e) {
            std::cerr << "[ERROR] -- IdentifierMapping::querySequenceIDConst -- sequence " + sequence
                         + " and/or genome " + genome + " not found" << std::endl;
            throw e;
        }
    }
    //! Returns the ID of a sequence (3)
    size_t querySequenceIDConst(SequenceGenomeTuple const & sequence) {
        try {
            return sequenceTupleToID_.at(sequence);
        } catch (std::exception const & e) {
            std::cerr << "[ERROR] -- IdentifierMapping::querySequenceIDConst -- sequence " + sequence.sequence
                         + " and/or genome " + std::to_string(sequence.gid) + " not found" << std::endl;
            throw e;
        }
    }
    //! Only for legacy code that is too cumbersome to update
    size_t querySequenceID_LEGACY(std::string const & sequence) const {
        for (auto&& gid : genomeIDToName_) {
            auto t = SequenceGenomeTuple{gid.first, sequence};
            if (sequenceTupleToID_.find(t) != sequenceTupleToID_.end()) {
                return sequenceTupleToID_.at(t);
            }
        }
        throw std::runtime_error("[ERROR] -- IdentifierMapping::querySequenceID_LEGACY -- sequence not found");
    }
    //! Returns the identifier string of a sequence
    /*! If the query is not present in the vector, it throws */
    std::string const & querySequenceName(size_t sid) const {
        try {
            return sequenceIDToTuple_.at(sid).sequence;   // throws if ID not present
        } catch (std::exception const & e) {
            std::cerr << "[ERROR] -- IdentifierMapping::querySequenceName -- sequence " + std::to_string(sid) + " not found" << std::endl;
            throw e;
        }
    }
    auto const & querySequenceTuple(size_t sid) const {
        try {
            return sequenceIDToTuple_.at(sid);  // throws if ID not present
        } catch (std::exception const & e) {
            std::cerr << "[ERROR] -- IdentifierMapping::querySequenceTuple -- sequence " + std::to_string(sid) + " not found" << std::endl;
            throw e;
        }
    }
    //! Const getter for sequenceNameToID
    auto const & sequenceIDToTuple() const{ return sequenceIDToTuple_; }
    //! return true if genome has an ID assigned
    bool sequenceKnown(std::string const & sequence, std::string const & genome) const {
        if (genomeKnown(genome)) {
            return sequenceTupleToID_.find(SequenceGenomeTuple{genomeNameToID_.at(genome), sequence}) != sequenceTupleToID_.end();
        } else {
            return false;
        }
    }
    //! Const getter for sequenceNameToID
    auto const & sequenceTupleToID() const{ return sequenceTupleToID_; }
    //! Query genomeID from sequenceID
    auto genomeIDFromSequenceID(size_t sequenceID) const {
        try {
            return sequenceIDToTuple_.at(sequenceID).gid;
        } catch (std::exception const & e) {
            std::cerr << "[ERROR] -- IdentifierMapping::genomeIDFromSequenceID -- sequence " + std::to_string(sequenceID) + " not found" << std::endl;
            throw e;
        }
    }
    //! Implements operator<< for an IdentifierMapping object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, IdentifierMapping const & im) {
        out << "Genome Name,  ID:" << std::endl << im.genomeNameToID_ << std::endl;
        out << "Sequence Tuple => ID:" << std::endl;
        for (auto&& elem : im.sequenceTupleToID_) {
            out << elem.first << " => " << elem.second << std::endl;
        }
        return out;
    }

private:
    //! Inserts a genome name into the mapping, returns false if name was already known
    bool insertNewGenome(std::string const & genomeName) {
        // only do stuff if genome is unknown
        if (genomeNameToID_.find(genomeName) == genomeNameToID_.end()) {
            genomeNameToID_.insert({genomeName, nextGenomeID_});
            genomeIDToName_.insert({nextGenomeID_, genomeName});
            ++nextGenomeID_;
            return true;
        }
        return false;
    }
    //! Inserts a sequence, returns false if sequence was already known
    bool insertNewSequence(std::string const & sequenceName, size_t gid) {
        if (genomeIDToName_.find(gid) == genomeIDToName_.end()) {
            throw std::runtime_error("[ERROR] -- IdentifierMapping::insertNewSequence -- genome ID not known");
        }
        auto sequenceTuple = SequenceGenomeTuple{gid, sequenceName};
        // only do stuff if sequence is unknown
        if (sequenceTupleToID_.find(sequenceTuple) == sequenceTupleToID_.end()) {
            sequenceTupleToID_.insert({sequenceTuple, nextSequenceID_});
            sequenceIDToTuple_.insert({nextSequenceID_, sequenceTuple});
            ++nextSequenceID_;
            return true;
        }
        return false;
    }

    tsl::hopscotch_map<size_t, std::string> genomeIDToName_;
    tsl::hopscotch_map<std::string, size_t> genomeNameToID_;
    size_t nextGenomeID_;
    size_t nextSequenceID_;
    tsl::hopscotch_map<size_t, SequenceGenomeTuple> sequenceIDToTuple_;
    tsl::hopscotch_map<SequenceGenomeTuple, size_t,
                       SequenceGenomeTupleHash> sequenceTupleToID_;
};

#endif // IDENTIFIERMAPPING_H
