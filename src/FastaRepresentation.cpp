#include "FastaRepresentation.h"

namespace fs = std::experimental::filesystem;

std::string FastaRepresentation::createRandomSequence(size_t length,
                                                      std::uniform_int_distribution<uint8_t> & unif,
                                                      std::default_random_engine & rng) {
    auto bases = std::array<char, 4>{'A', 'C', 'G', 'T'};
    auto seqVector = std::vector<char>(length);
    for (size_t i = 0; i < length; ++i) {
        seqVector[i] = bases.at(unif(rng)); // randomly concatenate bases
    }
    auto seq = std::string(seqVector.begin(), seqVector.end());
    return seq;
}



void FastaRepresentation::readFile(std::string const & fastaFile) {
    fs::path path(fastaFile);
    // Check if input file exists and is valid
    if (!fs::is_regular_file(path)) {
        auto msg = "[ERROR] -- FastaRepresentation -- " + fastaFile + " is not a regular file";
        throw std::runtime_error(msg);
    }
    // Try to open file (non binary mode)
    std::ifstream inputStream(path);
    if (!inputStream.is_open()) {
        auto msg = "[ERROR] -- FastaRepresentation -- Failed to open " + fastaFile;
        throw std::runtime_error(msg);
    }
    // If file is opened successfully, try to read contents
    std::string line;
    std::string currentHead;
    size_t lineCount{0};
    // getline returns inputStream and as https://en.cppreference.com/w/cpp/io/basic_ios/operator_bool,
    //     body only gets executed if getline was successful (i.e. no fail or read after eof)
    while ( std::getline(inputStream, line) ) {
        ++lineCount;
        chomp(line);    // strip newline sequences from line
        if (line.size() > 0) {  // skip empty lines
            if (line.at(0) == '>') { // new sequence
                line.erase(line.begin());   // remove ">"
                currentHead = line;
                // new (empty) sequence entry
                headToSeq_.insert({currentHead, FastaSequence("", currentHead, filename_)});
            } else if (line.at(0) == ';') { // skip comment line
                continue;
            } else {    // sequence line
                if (currentHead.empty()) {
                    auto msg = "[ERROR] -- FastaRepresentation -- Encountered non-header, non-comment line"
                               "that does not belong to any sequence header in '" + fastaFile +
                               "' in line " + std::to_string(lineCount);
                    throw std::runtime_error(msg);
                }
                headToSeq_.at(currentHead).append(line);    // create contiguous sequence from respective entries
            }
        }
    }
    if ( !inputStream.eof() && inputStream.fail() ) {
        auto msg = "[ERROR] -- FastaRepresentation -- Failed to read " + fastaFile;
        throw std::runtime_error(msg);
    }
}



FastaRepresentation::FastaRepresentation(std::string const & fastaFile)
    : artificialHeads_{}, filename_(genomeFromFilename(fastaFile)), headToSeq_{} {
    readFile(fastaFile);
}



FastaRepresentation::FastaRepresentation(std::string const & fastaFile, size_t artificialSequenceLength)
    : artificialHeads_{}, filename_{genomeFromFilename(fastaFile)}, headToSeq_{} {
    readFile(fastaFile);

    // create artificial sequence
    auto head = artificialHeader(fastaFile, artificialSequenceLength);
    std::random_device rd;
    std::default_random_engine rng(rd());
    std::uniform_int_distribution<uint8_t> unif(0, 3);
    auto seq = createRandomSequence(artificialSequenceLength, unif, rng);
    headToSeq_.insert({head, FastaSequence(seq, head, filename_)});
    artificialHeads_.insert(head);
}



FastaRepresentation::FastaRepresentation(std::string const & fastaFile, size_t artificialSizeFactor, DynamicallyGenerateArtificialSequences)
    : artificialHeads_{}, filename_{genomeFromFilename(fastaFile)}, headToSeq_{} {
    readFile(fastaFile);

    std::random_device rd;
    std::default_random_engine rng(rd());
    std::uniform_int_distribution<uint8_t> unif(0, 3);
    FastaMapType artificialMap;
    size_t count = 0;
    for (auto&& elem : headToSeq_) {
        for (size_t i = 0; i < artificialSizeFactor; ++i) {
            auto length = elem.second.sequence().size();
            auto head = artificialHeader(fastaFile, length, count);
            auto seq = createRandomSequence(length, unif, rng);
            artificialMap.insert({head, FastaSequence(seq, head, filename_)});
            artificialHeads_.insert(head);
            ++count;
        }
    }
    headToSeq_.insert(artificialMap.begin(), artificialMap.end());
}



void FastaRepresentation::chomp(std::string& line) {
    if (line.size() > 0) {
        // chomp \n and \r if present ( https://en.wikipedia.org/wiki/Newline#Representation )
        while (line.back() == '\n' || line.back() == '\r') {
            line.pop_back();
        }
    }
}



void FastaRepresentation::writeArtificialSequences(std::ofstream & os) const {
    if (!os.good()) { throw std::runtime_error("[ERROR] -- FastaRepresentation::writeArtificialSequences -- Cannot write to file"); }
    for (auto&& head : artificialHeads_) {
        auto& fa = headToSeq_.at(head);
        os << fa << std::endl;
    }
}
