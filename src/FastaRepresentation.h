#ifndef FASTAREPRESENTATION_H
#define FASTAREPRESENTATION_H

#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>

#include "Configuration.h"
#include "IdentifierMapping.h"
#include "StrongType.h"

namespace fs = std::filesystem;

//! Holds a single sequence from a fasta file, use with FastaRepresentation
class FastaSequence {
public:
    //! c'tor
    /*! \param sequence Sequence
     * \param sequenceName Name (header) of the sequence
     * \param genomeName Name of the genome from which the sequence is taken
     *
     * \details Initializes the members */
    FastaSequence(std::string const & sequence,
                  std::string const & sequenceName,
                  std::string const & genomeName)
        : genomeName_{genomeName},
          sequence_{sequence},
          sequenceName_{sequenceName} {}
    //! Append something to the sequence string
    void append(std::string const & extension) {
        sequence_.append(extension);
    }
    //! Getter for member \c genomeName_
    auto const & genomeName() const { return genomeName_; }
    bool operator==(FastaSequence const & rhs) const {
        return genomeName_ == rhs.genomeName_
                && sequence_ == rhs.sequence_
                && sequenceName_ == rhs.sequenceName_;
    }
    //! Getter for member \c sequence_
    auto const & sequence() const { return sequence_; }
    //! Getter for member \c sequenceName_
    auto const & sequenceName() const { return sequenceName_; }
    //! Implements operator<< for an FastaSequence object for use with \c std::ofstream (use for file output)
    friend std::ofstream & operator<<(std::ofstream & out, FastaSequence const & fasta) {
        out << ">" << fasta.sequenceName_ << std::endl;
        out << fasta.sequence_ << std::endl;
        return out;
    }
    //! Implements operator<< for an FastaSequence object for use with \c std::ostream (use for std::cout)
    friend std::ostream & operator<<(std::ostream & out, FastaSequence const & fasta) {
        out << ">" << fasta.sequenceName_ << std::endl;
        if (fasta.sequence_.size() > 30) {
            out << fasta.sequence_.substr(0,30) << "..." << std::endl;
        } else {
            out << fasta.sequence_ << std::endl;
        }
        return out;
    }

private:
    //! Stores the genome name
    std::string genomeName_;
    //! Stores the sequence
    std::string sequence_;
    //! Stores the sequence name
    std::string sequenceName_;
};



// Strong Types to distinguish constructor calls
using FastaGenomeName = NamedType<std::string, struct FastaGenomenameTag>;
using FastaFileName = NamedType<std::string, struct FastaFileNameTag>;

//! Holds contents of a fasta file, i.e. a mapping from fasta headers to the respective sequences
class FastaRepresentation {
public:
    using FastaMapType = std::unordered_map<std::string, FastaSequence>;
    static struct DynamicallyGenerateArtificialSequences{} dynamicallyGenerateArtificialSequences;

    // Factory functions
    //! Return a genome name from a file path
    static auto genomeFromFilename(std::string const & fastaFile) {
        fs::path path(fastaFile);
        auto filename = path.filename();
        filename.replace_extension("");
        return filename.string();
    }
    //! Creates a header for an artificial sequence
    static auto artificialHeader(std::string const & fastaFile, size_t artificialSequenceLength, size_t count = 0) {
        auto header = genomeFromFilename(fastaFile) + "|" + std::to_string(artificialSequenceLength) + "|artificial|" + std::to_string(count);
        return header;
    }

    //! c'tor (1)
    /*! \param genomeName Name used as genome name for this sequence collection
     *
     * \details Creates an empty FastaRepresentation under the given genome name */
    FastaRepresentation(FastaGenomeName const & genomeName)
        : artificialHeads_{}, filename_{genomeName.get()}, headToSeq_{} {}
    //! c'tor (2)
    /*! \param fastaFile Fasta file to read
     *
     * \details Reads the file from disk and stores its contents */
    FastaRepresentation(FastaFileName const & fastaFile);
    //! c'tor (3)
    /*! \param fastaFile Fasta file to read
     * \param artificialSequenceLength Length of the artificial sequence
     *
     * \details Reads the file from disk and stores its contents and creates
     * an additional random ACGT-sequence of length \c artificialSequenceLength */
    FastaRepresentation(FastaFileName const & fastaFile, size_t artificialSequenceLength);
    //! c'tor (4)
    /*! \param fastaFile Fasta file to read
     * \param DynamicallyGenerateArtificialSequences Tag to signal that artificial sequences are
     * to be generated dynamically to match the lenghts of the real sequences from \c fasta
     *
     * \details Reads the file from disk and stores its contents and creates
     * additional random ACGT-sequences that match the number and respective lengths of the
     * real sequences in the file */
    FastaRepresentation(FastaFileName const & fastaFile, size_t artificialSizeFactor, DynamicallyGenerateArtificialSequences);
    //! Add a new sequence to this FastaRepresentation
    void addSequence(std::string sequenceName, std::string const & sequence) {
        if (sequenceName.at(0) == '>') { sequenceName.erase(sequenceName.begin()); }    // remove ">"
        headToSeq_.insert({sequenceName, FastaSequence(sequence, sequenceName, filename_)});
    }
    //! Returns the FastaSequence mapped to \c header, throws if \c header is unknown
    auto const & fastaSequence(std::string const & header) const { return headToSeq_.at(header); }
    //! Getter for member \c filename_
    auto const & genome() const { return filename_; }
    //! Return a vector of sequence headers present in the fasta file
    auto headers() const {
        std::vector<std::string> headers;
        for (auto&& elem : headToSeq_) { headers.emplace_back(elem.first); }
        return headers;
    }
    //! Const reference to \c headToSeq_ private member
    auto const & headerToSequence() const { return headToSeq_; }
    //! Number of sequences in fasta file
    auto numSequences() const { return headToSeq_.size(); }
    //! Returns the sequence mapped to \c header, throws if \c header is unknown
    auto const & sequence(std::string const & header) const { return headToSeq_.at(header).sequence(); }
    //! Return number of sequences
    auto size() const { return headToSeq_.size(); }
    //! Write artificial sequences to file
    void writeArtificialSequences(std::ofstream & os) const;

private:
    //! Factory function to create a random sequence
    std::string createRandomSequence(size_t length,
                                     std::uniform_int_distribution<uint8_t> & unif,
                                     std::default_random_engine & rng);
    //! Factory function to read the file content from disc
    void readFile(std::string const & fastaFile);

    //! Remove newline and carriage return at the end of \c line in place
    void chomp(std::string& line);

    //! Store heads of artificial sequences
    std::set<std::string> artificialHeads_;
    //! Store the filename (without path and ending if any) of the input file
    std::string filename_;
    //! Stores a mapping from a fasta header to the respective sequence
    FastaMapType headToSeq_;
};

#endif // FASTAREPRESENTATION_H
