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
    Kmer(std::size_t K) : data(0), K(K) {}
    Kmer(const std::string &kmer) {
        data = 0;
        for (std::size_t i = 0; i < kmer.size(); i++) {
            if (i != 0) {
                data <<= 2;
            }
            data |= char_to_nucleotide(kmer[i]);
        }
        K = kmer.size();
    }
    void roll(char c) {
        data <<= 2;
        data |= char_to_nucleotide(c);
        data &= ((1ULL << (2 * K)) - 1);
    }
    /**
     * @brief Get the nucleotide at position i (0-indexed, from right to left)
     * @param i The position of the nucleotide to get
     * @return The nucleotide at position i
     */
    Nucleotide get(std::size_t i) const {
        if (i >= K) {
            return Nucleotide::N;
        }
        return static_cast<Nucleotide>((data >> (i * 2)) & 0b11);
    }
    Nucleotide last() const { return get(K - 1); }
    std::size_t size() const { return K; }
    void reset() { data = 0; }
    bool operator==(const Kmer &other) const {
        return data == other.data && K == other.K;
    }

  private:
    std::size_t K;
    std::uint64_t data;
};

#endif
