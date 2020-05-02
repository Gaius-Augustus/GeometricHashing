#ifndef REVERSECOMPLEMENT_H
#define REVERSECOMPLEMENT_H

#include <string>

inline char complement (char base) {
    switch(base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}



inline std::string reverseComplement (std::string const & fwd) {
    std::string revComp(fwd);
    auto last = fwd.length() - 1;
    for (size_t i = 0; i < fwd.length(); i++) {
        revComp[i] = complement(fwd[(last - i)]);
    }

    return revComp;
}

#endif // REVERSECOMPLEMENT_H
