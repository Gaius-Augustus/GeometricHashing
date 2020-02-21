#ifndef OPTIMALSPACEDSEEDS_H
#define OPTIMALSPACEDSEEDS_H

#include <exception>
#include <string>

#include "SpacedSeedMask.h"

inline std::vector<SpacedSeedMask> optimalSpacedSeeds(size_t weight) {
    switch (weight) {
    case 3:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("1101")};
        break;
    case 4:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("1100101")};
        break;
    case 5:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("10010111")};
        break;
    case 6:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("1100101011")};
        break;
    case 7:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11101010011")};
        break;
    case 8:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11001001010111")};
        break;
    case 9:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111000101011011")};
        break;
    case 10:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("1101100011010111")};
        break;
    case 11:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111001011001010111")};
        break;
    case 12:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111011001011010111")};
        break;
    case 13:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11101110010110010111")};
        break;
    case 14:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111101001101001110111")};
        break;
    case 15:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11110010101011001101111")};
        break;
    case 16:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111100110101011001101111")};
        break;
    case 17:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("1111010110011001101011111")};
        break;
    case 18:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111101100101101010110011111")};
        break;
    case 19:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11111011001110101101011111")};
        break;
    case 20:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11110110100111010011101110111")};
        break;
    case 21:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111110101100111001011101101111")};
        break;
    case 22:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111110111001110101011011011111")};
        break;
    case 23:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111111001101011010111001101011111")};
        break;
    case 24:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("1111110101011100111001011011011111")};
        break;
    case 25:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111110111011100111010110110111111")};
        break;
    case 26:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11111101011100111011011010111011111")};
        break;
    case 27:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111111011010111010111001110110111111")};
        break;
    case 28:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("1111110110101110110011110101110111111")};
        break;
    case 29:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11111011101101110011110110101110111111")};
        break;
    case 30:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("11111011110110111010111011110110111111")};
        break;
    case 31:
        return std::vector<SpacedSeedMask>{SpacedSeedMask("111111011110111010111011011110110111111")};
        break;
    default:
        break;
    }
    throw std::runtime_error("[ERROR] -- optimalSpacedSeeds -- no pre-computed optimal seed for weight " + std::to_string(weight));
}



inline std::vector<SpacedSeedMask> optimalSpacedSeeds(size_t weight, size_t m) {
    if (m == 1) { return optimalSpacedSeeds(weight); }
    if (m == 2) {
        switch (weight) {
        case 3:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("101001"), SpacedSeedMask("1000000000011") };
            break;
        case 4:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1100101"), SpacedSeedMask("11000000101") };
            break;
        case 5:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("10110011"), SpacedSeedMask("101001000011") };
            break;
        case 6:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1100101011"), SpacedSeedMask("110001000001101") };
            break;
        case 7:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11101001101"), SpacedSeedMask("111000100001011") };
            break;
        case 8:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("110110100111"), SpacedSeedMask("10110001000100111") };
            break;
        case 9:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110010111"), SpacedSeedMask("110010100000010100111") };
            break;
        case 10:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111010010110111"), SpacedSeedMask("11101001000100010111") };
            break;
        case 11:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110100110111"), SpacedSeedMask("111010001100010010111") };
            break;
        case 12:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111010110100110111"), SpacedSeedMask("11110010001000101001111") };
            break;
        case 13:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1111011001011010111"), SpacedSeedMask("1110100100010011000101111") };
            break;
        case 14:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111010110110011001111"), SpacedSeedMask("1110110000101000100101001111") };
            break;
        case 15:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111101100110101011111"), SpacedSeedMask("1110110010100011000010110111") };
            break;
        case 16:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1111011011100101110111"), SpacedSeedMask("111101010001011000010011001111") };
            break;
        case 17:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111101011011001101011111"), SpacedSeedMask("111110101000100001010010011001111") };
            break;
        case 18:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1111101011010111001101111"), SpacedSeedMask("1111101001100010001010000110110111") };
            break;
        case 19:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11111010110111001011101111"), SpacedSeedMask("11110111000100100010100010101101111") };
            break;
        case 20:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111011110110110011101011111"), SpacedSeedMask("11111011000101001000110001011011111") };
            break;
        case 21:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11110110101011011001111011111"), SpacedSeedMask("11111011001001100011000010010101011111") };
            break;
        case 22:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11111101101101011100111011111"), SpacedSeedMask("111110101100010010100010100100111011111") };
            break;
        case 23:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111110101111011100110110111111"), SpacedSeedMask("11111011010001010100100111000110111111") };
            break;
        case 24:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11111011101101011100111011011111"), SpacedSeedMask("111110111000101100101001011000110111111") };
            break;
        }
    }
    if (m == 4) {
        switch (weight) {
        case 3:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("10011"), SpacedSeedMask("101000001"), SpacedSeedMask("10000000011"), SpacedSeedMask("1000010000001") };
            break;
        case 4:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("101011"), SpacedSeedMask("11000101"), SpacedSeedMask("1010000011"), SpacedSeedMask("10010000000011") };
            break;
        case 5:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11101001"), SpacedSeedMask("1001100011"), SpacedSeedMask("1010001011"), SpacedSeedMask("11000000100101") };
            break;
        case 6:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111010011"), SpacedSeedMask("110001001011"), SpacedSeedMask("11001000001011"), SpacedSeedMask("101010000001000011") };
            break;
        case 7:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1110101011"), SpacedSeedMask("110100100111"), SpacedSeedMask("1100010001001011"), SpacedSeedMask("110010000001010011") };
            break;
        case 8:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("110101011011"), SpacedSeedMask("111001000101011"), SpacedSeedMask("111000101000010011"), SpacedSeedMask("110100100000010000111") };
            break;
        case 9:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110100111"), SpacedSeedMask("11010000110010111"), SpacedSeedMask("11100010010000101011"), SpacedSeedMask("11010001000001000100111") };
            break;
        case 10:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11101010110111"), SpacedSeedMask("110110010100010111"), SpacedSeedMask("11101001000100001010011"), SpacedSeedMask("11010010000100001001000111") };
            break;
        case 11:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111011011010111"), SpacedSeedMask("111010100011001111"), SpacedSeedMask("1110100010010100100111"), SpacedSeedMask("11110001001000000100010111") };
            break;
        case 12:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110110101111"), SpacedSeedMask("1111010100011010111"), SpacedSeedMask("111100010010000100110111"), SpacedSeedMask("11101000100010000010010001111") };
            break;
        case 13:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1110110101001101111"), SpacedSeedMask("111010100010110000111011"), SpacedSeedMask("111010011100001001010111"), SpacedSeedMask("111100100010000001000100101111") };
            break;
        case 14:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111101101011101111"), SpacedSeedMask("11101100110010010101111"), SpacedSeedMask("111100011001010001001001111"), SpacedSeedMask("1110100101100000001000110010111") };
            break;
        case 15:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11101110110110101111"), SpacedSeedMask("1111101001010001100110111"), SpacedSeedMask("111101000100010010010100011111"), SpacedSeedMask("111100011001000100001000101001111") };
            break;
        case 16:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111110110011110101111"), SpacedSeedMask("1111001001010001100101011111"), SpacedSeedMask("111101010100001101001001100111"), SpacedSeedMask("11101100010010001000010100001101111") };
            break;
        case 17:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1111010111011101101111"), SpacedSeedMask("11110110001101001010001011111"), SpacedSeedMask("111101000101001001000101100011111"), SpacedSeedMask("11110010110000010001000100100001110111") };
            break;
        case 18:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11101111010110110111111"), SpacedSeedMask("11110110011000110101010011111"), SpacedSeedMask("111110010010101001001001100011111"), SpacedSeedMask("111110001100100100000101000100001011111") };
            break;
        case 19:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11110110111010111001101111"), SpacedSeedMask("1111001010011100101010011011111"), SpacedSeedMask("11111001100100001001101001010101111"), SpacedSeedMask("1111010100100010100010000011001001011111") };
            break;
        case 20:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11111011010111011101101111"), SpacedSeedMask("111110110101001100001100101011111"), SpacedSeedMask("111101100010100001100010010010101011111"), SpacedSeedMask("1111100101001010001000000101100001011001111") };
            break;
        case 21:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111110110111101011101101111"), SpacedSeedMask("111110110001100101011010100111111"), SpacedSeedMask("11111010100011010010000110010011011111"), SpacedSeedMask("1111100110100010000010100001100101001011111") };
            break;
        case 22:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111110111101101110111011111"), SpacedSeedMask("111111011000111000100101001101011111"), SpacedSeedMask("11101011101100001011010001011001101111"), SpacedSeedMask("11111100101000100100100000110010100011011111") };
            break;
        case 23:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11111011110111011101101101111"), SpacedSeedMask("11111011011000101001110011010111111"), SpacedSeedMask("11101101100001111101001000110101011111"), SpacedSeedMask("1111011100011010100000011001001001011011111") };
            break;
        case 24:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111110111101011111011011110111"), SpacedSeedMask("1111101101011010110001001110010111111"), SpacedSeedMask("111101110011001000101110100101110101111"), SpacedSeedMask("111111001100100101000100000110010101001111111") };
            break;
        }
    }
    if (m == 8) {
        switch (weight) {
        case 3:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1011"), SpacedSeedMask("1010000001"), SpacedSeedMask("1000000000011"), SpacedSeedMask("10010000000001"), SpacedSeedMask("100010000000001"), SpacedSeedMask("100001000000001"), SpacedSeedMask("100000100000001"), SpacedSeedMask("100000001000001") };
            break;
        case 4:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("10111"), SpacedSeedMask("10001000011"), SpacedSeedMask("10010000000101"), SpacedSeedMask("100000000100101"), SpacedSeedMask("1000001000000101"), SpacedSeedMask("1001000010000001"), SpacedSeedMask("1000100000000011"), SpacedSeedMask("1000000010000011") };
            break;
        case 5:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("110111"), SpacedSeedMask("110001011"), SpacedSeedMask("101001000011"), SpacedSeedMask("1010001001001"), SpacedSeedMask("10100100000011"), SpacedSeedMask("100100010000101"), SpacedSeedMask("110000100010001"), SpacedSeedMask("1010000010000000011") };
            break;
        case 6:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11010111"), SpacedSeedMask("11100100101"), SpacedSeedMask("110010001011"), SpacedSeedMask("1011000010011"), SpacedSeedMask("11000010100011"), SpacedSeedMask("1010000100010011"), SpacedSeedMask("1100100000000010101"), SpacedSeedMask("11000010000000001101") };
            break;
        case 7:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("110101111"), SpacedSeedMask("110101000111"), SpacedSeedMask("11100100010011"), SpacedSeedMask("1100101000010101"), SpacedSeedMask("11010001000001011"), SpacedSeedMask("11000100100000100011"), SpacedSeedMask("110100000000100001011"), SpacedSeedMask("1010010000010000010011") };
            break;
        case 8:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("110100111011"), SpacedSeedMask("11101010001011"), SpacedSeedMask("110101001000111"), SpacedSeedMask("110010010010000111"), SpacedSeedMask("1101000001000101011"), SpacedSeedMask("1100101000010010011"), SpacedSeedMask("110101000010000010011"), SpacedSeedMask("11100010000000100001011") };
            break;
        case 9:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11110110111"), SpacedSeedMask("11101010001111"), SpacedSeedMask("1101100010010111"), SpacedSeedMask("11100001010001010011"), SpacedSeedMask("110010010000100011011"), SpacedSeedMask("101011000001001000111"), SpacedSeedMask("1100101000100001001011"), SpacedSeedMask("111000010010000001000111") };
            break;
        case 10:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111110110111"), SpacedSeedMask("11110001010011011"), SpacedSeedMask("110101100001010111"), SpacedSeedMask("11010100010011000111"), SpacedSeedMask("1101000010100101000111"), SpacedSeedMask("11100100100000100101011"), SpacedSeedMask("11100010001000000100100111"), SpacedSeedMask("11101000000100010001001011") };
            break;
        case 11:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11101011011111"), SpacedSeedMask("11101001101010111"), SpacedSeedMask("111001010110000010111"), SpacedSeedMask("111010100000110100111"), SpacedSeedMask("1101100010010001101011"), SpacedSeedMask("111000101001001000100111"), SpacedSeedMask("110110001100010001001011"), SpacedSeedMask("111100010000010000100010111") };
            break;
        case 12:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111011101101111"), SpacedSeedMask("11101100001100110111"), SpacedSeedMask("111100011010001010111"), SpacedSeedMask("111010010010100100001111"), SpacedSeedMask("111001010101000011001011"), SpacedSeedMask("11101000110000010010100111"), SpacedSeedMask("1101101000001000100001010111"), SpacedSeedMask("11100010010000100010010010111") };
            break;
        case 13:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1111011011101111"), SpacedSeedMask("11101011100010110111"), SpacedSeedMask("111101000101101001111"), SpacedSeedMask("1101100011010000101010111"), SpacedSeedMask("110101110000100010000110111"), SpacedSeedMask("111010010010001010001001111"), SpacedSeedMask("1110110001000011000100100111"), SpacedSeedMask("111001010010010000100100010111") };
            break;
        case 14:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111101101011101111"), SpacedSeedMask("111101011100010110111"), SpacedSeedMask("111010100001011001101111"), SpacedSeedMask("1110110001101000101010111"), SpacedSeedMask("11101001100011001100011011"), SpacedSeedMask("111001011001000100101001111"), SpacedSeedMask("11101010010001000001000100101111"), SpacedSeedMask("11101101000010010000010000111011") };
            break;
        case 15:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111111011101101111"), SpacedSeedMask("1111101001011100110111"), SpacedSeedMask("1101010111000100011101111"), SpacedSeedMask("11110011000110100011010111"), SpacedSeedMask("1111001001000101101000101111"), SpacedSeedMask("11101101001010001010010010111"), SpacedSeedMask("11110100010100001000100110001111"), SpacedSeedMask("111101010000100100000010010110111") };
            break;
        case 16:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1111110111101110111"), SpacedSeedMask("111011100011100101101111"), SpacedSeedMask("1111010011000101110110111"), SpacedSeedMask("1111010010011010010001011111"), SpacedSeedMask("11110011001010000101100001011011"), SpacedSeedMask("110111001000010100100010110001111"), SpacedSeedMask("11110010101000100000001100101001111"), SpacedSeedMask("11110110010010000010100001001100111") };
            break;
        case 17:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11111101111011101111"), SpacedSeedMask("1111101011001011100101111"), SpacedSeedMask("111010110001101001101011111"), SpacedSeedMask("111101001101000010101001110111"), SpacedSeedMask("11011101000010110010000110101111"), SpacedSeedMask("1111100100101000100011010001001111"), SpacedSeedMask("111100110010000100101000110010001111"), SpacedSeedMask("111100101000110000100010000110110111") };
            break;
        case 18:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1111101110110111101111"), SpacedSeedMask("111110110101100101100110111"), SpacedSeedMask("1110110100111010000110001011111"), SpacedSeedMask("111101010001100001010011110010111"), SpacedSeedMask("111101010100000110110001001101111"), SpacedSeedMask("1110011011001010000001010101001011011"), SpacedSeedMask("1111010000110010100100100010010101111"), SpacedSeedMask("1111000110001001000110010100110001111") };
            break;
        case 19:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111110110110111011110111"), SpacedSeedMask("11111001110010101100101011111"), SpacedSeedMask("111011101010011000011110110111"), SpacedSeedMask("111010111001000101110010010011111"), SpacedSeedMask("11110100100110100010001010110011111"), SpacedSeedMask("111101011000100101000101000010001110111"), SpacedSeedMask("111100100111000000101001001101000110111"), SpacedSeedMask("111110100000101001010000101001011001111") };
            break;
        case 20:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111101111011011111011111"), SpacedSeedMask("11110110010111011001010111111"), SpacedSeedMask("111110011101010010111001101111"), SpacedSeedMask("11111100100010110001110100100101111"), SpacedSeedMask("11110101011100100100000110011011111"), SpacedSeedMask("1111010110001100000011001001010010011111"), SpacedSeedMask("1111100100101001100010000110001010101111"), SpacedSeedMask("1110110011010000010100110100001010110111") };
            break;
        case 21:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11111011101110110111101111"), SpacedSeedMask("1111101001101001100111010111111"), SpacedSeedMask("1111011001010101101011100100011111"), SpacedSeedMask("111010101100111110000010100110110111"), SpacedSeedMask("111110111001010000010100100011011010111"), SpacedSeedMask("111101101001100011001000000011100010101111"), SpacedSeedMask("111110001011010000001100010100001101110111"), SpacedSeedMask("111101001100000101100001101010001001011111") };
            break;
        case 22:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("111011011101111101111011111"), SpacedSeedMask("11111001101111001001101010111111"), SpacedSeedMask("1111101101000100111010110110011111"), SpacedSeedMask("1111001101010111000110001001011011111"), SpacedSeedMask("1110111101100000010101100101010010110111"), SpacedSeedMask("1111100101010110001001000010100011000111111"), SpacedSeedMask("1111001011000010110100010000111000101110111"), SpacedSeedMask("1110111010100001100000011011000110010110111") };
            break;
        case 23:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("1111011110111101110110111111"), SpacedSeedMask("11111011001101010111110011101111"), SpacedSeedMask("11110100111100111110010100101011111"), SpacedSeedMask("111011110000111011000101011101011111"), SpacedSeedMask("11111101010100100001101110010011011111"), SpacedSeedMask("111110100110011010100000100010110100111111"), SpacedSeedMask("111110101100101001000010011000001100111001111"), SpacedSeedMask("111111010011000010001011000011000100101011111") };
            break;
        case 24:
            return std::vector<SpacedSeedMask>{ SpacedSeedMask("11111101101110111011110111111"), SpacedSeedMask("1111100111010110100101110110111111"), SpacedSeedMask("11110110101101101111100100011101111"), SpacedSeedMask("11111001111001000110010110101101011111"), SpacedSeedMask("11110111001100001111010001010010010111111"), SpacedSeedMask("11110100110001110010101011100001110101111"), SpacedSeedMask("11111101010010011000010010001011100011101111"), SpacedSeedMask("1111101001101010000001100011000100110110011111") };
            break;
        }
    }

    throw std::runtime_error("[ERROR] -- optimalSpacedSeeds -- no pre-computed optimal seed for weight " + std::to_string(weight) + " and seed set size " + std::to_string(m));
}
#endif // OPTIMALSPACEDSEEDS_H
