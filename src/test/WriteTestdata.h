#ifndef WRITETESTDATA_H
#define WRITETESTDATA_H

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

class WriteTestdata {
public:
    WriteTestdata(std::string const & content)
        : hashfun_{}, filename_{"tmp_testFile_" + std::to_string(hashfun_(content))} {

        fs::path tempFile(filename_);
        if (fs::exists(tempFile)) {
            std::cerr << "[WARNING] -- WriteTestdata -- File '" << filename_ << "' already exists, deleting... " << std::endl;
            deleteFile();
        }
        if (fs::exists(tempFile)) { throw std::runtime_error("[ERROR] -- WriteTestdata -- File already exists"); }

        std::ofstream os(tempFile);
        if (os.good()) {
            os << content;
            os.close();
        } else {
            throw std::runtime_error("[ERROR] -- WriteTestdata -- Stream corrupt");
        }
    }
    void deleteFile() const {
        fs::path tempFile(filename_);
        if (fs::exists(tempFile)) {
            fs::remove(tempFile);
        } else {
            throw std::runtime_error("[ERROR] -- WriteTestdata::deleteFile -- File not present");
        }
    }
    auto const & filename() const { return filename_; }
private:
    std::hash<std::string> hashfun_;
    std::string const filename_;
};

#endif // WRITETESTDATA_H
