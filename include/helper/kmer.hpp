#ifndef KMER_HPP
#define KMER_HPP

#include <cstdint>
#include <stdexcept>

enum Nucleotide { A = 0, C = 1, G = 2, T = 3, N = 4 };

constexpr inline Nucleotide char_to_nucleotide(char c) {
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

class Kmer {
  public:
    using data_t = std::uint64_t;
    explicit Kmer(std::size_t K) : _data(0), K(K), n_count(0) {}
    explicit Kmer(const std::string &kmer) {
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
    void roll(char c) {
        _data <<= 2;
        _data |= char_to_nucleotide(c);
        _data &= ((1ULL << (2 * K)) - 1);
        n_count++;
    }
    /**
     * @brief Get the nucleotide at position i (0-indexed, from right to left)
     * @param i The position of the nucleotide to get
     * @return The nucleotide at position i
     */
    Nucleotide get(std::size_t i) const {
        if (i >= K || i >= n_count) {
            return Nucleotide::N;
        }
        return static_cast<Nucleotide>((_data >> (i * 2)) & 0b11);
    }
    Nucleotide last() const { return get(K - 1); }
    std::size_t size() const { return K; }
    std::size_t available() const { return std::min(K, n_count); }
    void reset() {
        _data = 0;
        n_count = 0;
    }
    bool operator==(const Kmer &other) const {
        return other.K == K && _data == other._data;
    }
    const data_t &data() const { return _data; }

  private:
    std::size_t K;
    data_t _data;
    std::size_t n_count;
};

#endif
