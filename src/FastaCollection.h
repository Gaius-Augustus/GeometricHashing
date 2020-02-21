#ifndef FASTACOLLECTION_H
#define FASTACOLLECTION_H

#include <string>
#include <unordered_map>

#include "Configuration.h"
#include "FastaRepresentation.h"
#include "IdentifierMapping.h"

class FastaCollection {
public:
    //! c'tor (1)
    /*! Creates an empty collection */
    FastaCollection() : collection_{} {}
    //! c'tor (2)
    /*! \param config Shared pointer to program configuration object
     *
     * \details Creates and stores a FastaRepresentation for each input file
     * specified in \c config */
    FastaCollection(std::shared_ptr<Configuration const> config)
        : collection_{} {
        for (auto&& file : config->inputFiles()) {
            auto genomeName = FastaRepresentation::genomeFromFilename(file);
            if (/*config->artificialSequenceLengths().size() == 0 && */!config->dynamicArtificialSequences()) {
                collection_.emplace(genomeName, file);
            } else /*if (config->dynamicArtificialSequences())*/ {
                collection_.emplace(genomeName, FastaRepresentation(file,
                                                                    config->artificialSequenceSizeFactor(),
                                                                    FastaRepresentation::dynamicallyGenerateArtificialSequences));
            }
        }
    }
    //! Getter for member \c collection_
    auto const & collection() const { return collection_; }
    //! Insert new FastaRepresentation
    void emplace(std::string const & genomeName, std::string const & filename) {
        collection_.emplace(genomeName, filename);
    }
    //! Getter for the FastaRepresentation of \c genomeName
    auto const & fastaRepresentation(std::string const & genomeName) const {
        return collection_.at(genomeName);
    }
    auto numSequences() const {
        size_t n = 0;
        for (auto&& fa : collection_) {
            n += fa.second.numSequences();
        }
        return n;
    }
    void populateIdentifierMappingFromFastaCollection(IdentifierMapping & idMap) const {
        for (auto&& fasta : collection_) {
            idMap.queryGenomeID(fasta.first);
            for (auto&& seq : fasta.second.headerToSequence()) {
                idMap.querySequenceID(seq.second.sequenceName(), seq.second.genomeName());
            }
        }
    }
    FastaSequence const & fastaSequence(uint32_t sequenceID, IdentifierMapping const & idMap) const {
        auto sequenceName = idMap.querySequenceName(sequenceID);
        auto genomeID = idMap.genomeIDFromSequenceID(sequenceID);
        auto genomeName = idMap.queryGenomeName(genomeID);
        return collection_.at(genomeName).fastaSequence(sequenceName);
    }
    auto sequenceLength(uint32_t sequenceID, IdentifierMapping const & idMap) const {
        auto& seq = fastaSequence(sequenceID, idMap);
        return seq.sequence().size();
    }
private:
    //! genomeName to FastaRepresentation
    std::unordered_map<std::string,
                       FastaRepresentation> collection_;
};

#endif // FASTACOLLECTION_H
