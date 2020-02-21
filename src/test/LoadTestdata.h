#ifndef LOADTESTDATA_H
#define LOADTESTDATA_H

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../FastaCollection.h"
#include "../FastaRepresentation.h"
#include "../IdentifierMapping.h"
#include "GeometricHashingTestdata.h"

/* template <template<typename, typename> class SeedMapType,
          typename TwoBitKmerDataType,
          typename TwoBitSeedDataType> */
class LoadTestdata {
public:
    static struct GeometricHashingData{} geometricHashingData; // tag dispatch

    LoadTestdata()
        : fastaCollection_{std::make_shared<FastaCollection>()},
          genome0_{"hg38_orthologs"},
          genome1_{"mm10_orthologs"},
          idMap_{std::make_shared<IdentifierMapping>(genome0_)},    // id 0
          inputFiles_{"./testdata/hetGla2_orthologs.fa",
                      "./testdata/hg38_orthologs.fa",
                      "./testdata/macFas5_orthologs.fa",
                      "./testdata/mm10_orthologs.fa"} {
        idMap_->queryGenomeID(genome1_);    // id 1
        for (auto&& file : inputFiles_) {
            fastaCollection_->emplace(FastaRepresentation::genomeFromFilename(file), file);
        }
        fastaCollection_->populateIdentifierMappingFromFastaCollection(*idMap_);
        std::cout << "Current working dir: " << std::filesystem::current_path() << std::endl;
    }
    LoadTestdata(GeometricHashingData)
        : fastaCollection_{std::make_shared<FastaCollection>()},
          genome0_{kTestSpecies[0]},
          genome1_{kTestSpecies[1]},
          idMap_{std::make_shared<IdentifierMapping>(genome0_)},    // id 0
          inputFiles_{"./testdata/geometricHashing/species1.fasta",
                      "./testdata/geometricHashing/species2.fasta",
                      "./testdata/geometricHashing/species3.fasta",
                      "./testdata/geometricHashing/species4.fasta",
                      "./testdata/geometricHashing/species5.fasta"} {
        idMap_->queryGenomeID(genome1_);    // id 1
        idMap_->queryGenomeID(kTestSpecies[2]);    // id 2
        idMap_->queryGenomeID(kTestSpecies[3]);    // id 3
        idMap_->queryGenomeID(kTestSpecies[4]);    // id 4
        idMap_->querySequenceID(kTestSequences[0], kTestSpecies[0]);    // id 0
        idMap_->querySequenceID(kTestSequences[1], kTestSpecies[1]);    // id 1
        idMap_->querySequenceID(kTestSequences[2], kTestSpecies[1]);    // id 2
        idMap_->querySequenceID(kTestSequences[3], kTestSpecies[2]);    // id 3
        idMap_->querySequenceID(kTestSequences[4], kTestSpecies[2]);    // id 4
        idMap_->querySequenceID(kTestSequences[5], kTestSpecies[3]);    // id 5
        idMap_->querySequenceID(kTestSequences[6], kTestSpecies[3]);    // id 6
        idMap_->querySequenceID(kTestSequences[7], kTestSpecies[3]);    // id 7
        idMap_->querySequenceID(kTestSequences[8], kTestSpecies[4]);    // id 8
        idMap_->querySequenceID(kTestSequences[9], kTestSpecies[4]);    // id 9

        for (auto&& file : inputFiles_) {
            fastaCollection_->emplace(FastaRepresentation::genomeFromFilename(file), file);
        }
        std::cout << "Current working dir: " << std::filesystem::current_path() << std::endl;
    }

    auto fastaCollection() const { return fastaCollection_; }
    auto const & genome0() const { return genome0_; }
    auto const & genome1() const { return genome1_; }
    auto idMap() const { return idMap_; }
    auto const & inputFiles() const { return inputFiles_; }
//    auto seedMap() const { return seedMap_; }

private:
    std::shared_ptr<FastaCollection> fastaCollection_;
    std::string genome0_;
    std::string genome1_;
    std::shared_ptr<IdentifierMapping> idMap_;
    std::vector<std::string> inputFiles_;
};

#endif // LOADTESTDATA_H
