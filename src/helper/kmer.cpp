#include "helper/kmer.hpp"

#include <stdexcept>

using data_t = Kmer::data_t;

Nucleotide char_to_nucleotide(char c) {
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
    K = kmer.size();
    _data = 0;
    _rev_data = 0;
    for (std::size_t i = 0; i < kmer.size(); i++) {
        _data |= char_to_nucleotide(kmer[i]) << (2 * (K - 1 - i));
        _rev_data |= COMPLEMENT[char_to_nucleotide(kmer[i])] << (2 * i);
    }
    n_count = kmer.size();
}

void Kmer::roll(char c) {
    _data <<= 2;
    _data |= char_to_nucleotide(c);
    _data &= ((1ULL << (2 * K)) - 1);

    _rev_data >>= 2;
    _rev_data |= (data_t)COMPLEMENT[char_to_nucleotide(c)] << (2 * (K - 1));

    n_count++;
}

Nucleotide Kmer::get(std::size_t i, KmerRepr representation) const {
    if (i >= K || i >= n_count) {
        return Nucleotide::N;
    }
    return static_cast<Nucleotide>((data(representation) >> (i * 2)) & 0b11);
}
