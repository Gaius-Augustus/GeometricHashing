#ifndef CONFIGURATIONGENERATOR_H
#define CONFIGURATIONGENERATOR_H

#include <memory>
#include <string>
#include <vector>

#include "Configuration.h"

// https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
// split strings at each space and insert parts into argvStr
inline void addParameterFromString(std::vector<std::string> & argvStr, std::string const & str) {
    std::istringstream iss(str);
    argvStr.insert(argvStr.end(),
                   std::istream_iterator<std::string>(iss),
                   std::istream_iterator<std::string>());
}



inline auto configurationPtrFromParameters(std::vector<std::string> const & argvStr) {
    // convert into C++-style argv
    char** argv = new char*[argvStr.size() + 1];
    for (size_t i = 0; i < argvStr.size(); ++i) {
        argv[i] = new char[argvStr.at(i).size() + 1];   // +1 for null termination
        argvStr.at(i).copy(argv[i], argvStr.at(i).size());
        argv[i][argvStr.at(i).size()] = 0;  // null terminating
    }
    argv[argvStr.size()] = nullptr;

    //for (size_t i = 0; i < argvStr.size(); ++i) { std::cout << argv[i] << std::endl; }

    auto config = std::make_shared<Configuration>(argvStr.size(), argv);

    // cleanup raw pointer mess
    for (size_t i = 0; i < argvStr.size() + 1; ++i) { delete[] argv[i]; }
    delete[] argv;

    return config;
}



inline auto configurationFromParameters(std::vector<std::string> const & argvStr) {
    auto config = configurationPtrFromParameters(argvStr);
    return *config;
}



inline auto generateConfiguration(std::string const & str) {
    std::vector<std::string> argvStr;
    addParameterFromString(argvStr, str);
    return configurationPtrFromParameters(argvStr);
}

#endif // CONFIGURATIONGENERATOR_H
