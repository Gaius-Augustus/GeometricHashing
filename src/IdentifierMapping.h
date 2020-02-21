#ifndef IDENTIFIERMAPPING_H
#define IDENTIFIERMAPPING_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "prettyprint/prettyprint.hpp"



//! Stores a mapping of genome and sequence IDs to the identifier strings
/*! The reference genome always gets genome ID 0. The other IDs are assigned
 * by incrementing the last ID each time a new string is queried */
class IdentifierMapping {
public:
    //! Constructor
    /*! \param referenceGenome Identifier string of the reference genome
     *
     * \details Assigns genome ID 0 to the reference genome */
    IdentifierMapping(std::string const & referenceGenome)
        : genomeName_{}, genomeNameToID_{}, referenceGenome_{referenceGenome},
          referenceGenomeIDQueried_{0}, sequenceIDToGenomeID_{},sequenceName_{},
          sequenceNameToID_{} {
        // reference genome always gets ID '0'
        genomeName_.emplace_back(referenceGenome_);
        genomeNameToID_[referenceGenome_] = 0;
    }

    //! Copies state of the source IdentifierMapping
    void completeUpdate(IdentifierMapping const & source) {
        *this = source;
    }
    //! Returns the total number of genomes in the input data
    size_t numGenomes() const { return genomeName_.size(); }
    //! Returns the total number of sequences in the input data
    size_t numSequences() const { return sequenceName_.size(); }
    //! Memberwise comparison
    bool operator==(IdentifierMapping const & rhs) const {
        return genomeName_ == rhs.genomeName_
                && genomeNameToID_ == rhs.genomeNameToID_
                && referenceGenome_ == rhs.referenceGenome_
                && referenceGenomeIDQueried_ == rhs.referenceGenomeIDQueried_
                && sequenceName_ == rhs.sequenceName_
                && sequenceNameToID_ == rhs.sequenceNameToID_;
    }
    //! Returns identifier string of reference genome
    std::string const & refGenome() const { return referenceGenome_; }
    //! Returns the ID of a genome (1)
    /*! If the query is not present in the map, it is inserted now
     * Increments reference query count if parameter is query identifier string */
    size_t queryGenomeID(std::string const & genome) {
        if (genome == referenceGenome_) { ++referenceGenomeIDQueried_; }
        if (genomeNameToID_.find(genome) == genomeNameToID_.end()) {
            genomeName_.emplace_back(genome);
            genomeNameToID_[genome] = (genomeName_.size() - 1);
        }
        return genomeNameToID_.at(genome);
    }
    //! Returns the ID of a genome (2) -- throws it genome not known
    size_t queryGenomeIDConst(std::string const & genome) const {
        return genomeNameToID_.at(genome);
    }
    //! Returns the identifier string of a genome
    /*! If the query is not present in the vector, it throws */
    std::string const & queryGenomeName(size_t const id) const {
        return genomeName_.at(id);   // throws if ID not present
    }
    //! Returns the ID of a sequence (1)
    size_t querySequenceID(std::string const & sequence) const {
        return sequenceNameToID_.at(sequence);
    }
    //! Returns the ID of a sequence (2)
    /*! If the query is not present in the map, it is inserted now */
    size_t querySequenceID(std::string const & sequence, std::string const & genome) {
        if (sequenceNameToID_.find(sequence) == sequenceNameToID_.end()) {
            sequenceName_.emplace_back(sequence);
            sequenceNameToID_[sequence] = (sequenceName_.size() - 1);
            auto genomeID = queryGenomeID(genome);
            auto sequenceID = querySequenceID(sequence);
            sequenceIDToGenomeID_.insert({sequenceID, genomeID});
        }
        return sequenceNameToID_.at(sequence);
    }
    //! Returns the identifier string of a sequence
    /*! If the query is not present in the vector, it throws */
    std::string const & querySequenceName(uint32_t const id) const {
        return sequenceName_.at(id);   // throws if ID not present
    }
    //! Getter for member \c sequenceIDToGenomeID_
    auto const & sequenceIDToGenomeID() const { return sequenceIDToGenomeID_; }
    //! Query genomeID from sequenceID
    auto genomeIDFromSequenceID(size_t sequenceID) const { return sequenceIDToGenomeID_.at(sequenceID); }
    //! Implements operator<< for an IdentifierMapping object for use with \c std::ostream
    friend std::ostream & operator<<(std::ostream & out, IdentifierMapping const & im) {
        out << "Reference: " << im.referenceGenome_ << std::endl;
        out << "Genome Name => ID:" << std::endl << im.genomeNameToID_ << std::endl;
        out << "Sequence Name => ID:" << std::endl << im.sequenceNameToID_ << std::endl;
        out << "Sequence ID => Genome ID:" << std::endl << im.sequenceIDToGenomeID_ << std::endl;
        return out;
    }

private:
    //! Stores genome identifier strings
    std::vector<std::string> genomeName_;
    //! Stores mapping from genome identifier strings to their IDs
    std::unordered_map<std::string, size_t> genomeNameToID_;
    //! Stores reference genome identifier string
    std::string referenceGenome_;
    //! Count of how often the reference genome ID was queried
    size_t referenceGenomeIDQueried_;
    //! Stores mapping from a sequence ID to its respective genome ID
    std::unordered_map<size_t, size_t> sequenceIDToGenomeID_;
    //! Stores sequence identifier strings
    std::vector<std::string> sequenceName_;
    //! Stores mapping from sequence identifier strings to their IDs
    std::unordered_map<std::string, size_t> sequenceNameToID_;
};

#endif // IDENTIFIERMAPPING_H
