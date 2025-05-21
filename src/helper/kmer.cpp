#include "helper/kmer.hpp"

#include <stdexcept>

using data_t = Kmer::data_t;

constexpr Nucleotide reverse_complement[] = {T, G, C, A};

constexpr Nucleotide char_to_nucleotide(char c) {
    switch (c) {
    case 'A':
    case 'a':
        return A;
    case 'C':
    case 'c':
        return C;
    case 'G':
    case 'g':
        return G;
    case 'T':
    case 't':
        return T;
    case 'N':
        return N;
    default:
        throw std::invalid_argument("Invalid nucleotide character");
    }
}

Kmer::Kmer(const std::string &kmer) {
    _data = 0;
    for (std::size_t i = 0; i < kmer.size(); i++) {
        if (i != 0) {
            _data <<= 2;
        }
        _data |= char_to_nucleotide(kmer[i]);
    }
    K = kmer.size();
    n_count = kmer.size();
}

void Kmer::roll(char c) {
    _data <<= 2;
    _data |= char_to_nucleotide(c);
    _data &= ((1ULL << (2 * K)) - 1);
    n_count++;
}

Nucleotide Kmer::get(std::size_t i) const {
    if (i >= K || i >= n_count) {
        return Nucleotide::N;
    }
    return static_cast<Nucleotide>((_data >> (i * 2)) & 0b11);
}

data_t Kmer::reverse_complement(data_t data) const {
    data = (data >> 2 & 0x3333333333333333ULL) |
           (data << 2 & 0xCCCCCCCCCCCCCCCCULL);
    data = (data >> 4 & 0x0F0F0F0F0F0F0F0FULL) |
           (data << 4 & 0xF0F0F0F0F0F0F0F0ULL);
    data = (data >> 8 & 0x00FF00FF00FF00FFULL) |
           (data << 8 & 0xFF00FF00FF00FF00ULL);
    data = (data >> 16 & 0x0000FFFF0000FFFFULL) |
           (data << 16 & 0xFFFF0000FFFF0000ULL);
    data = (data >> 32) | (data << 32);

    data_t mask = (1ULL << (2 * K)) - 1;
    data >>= 8 * sizeof(data) - 2 * K;

    return ~data & mask;
}
