#ifndef KMER_HPP
#define KMER_HPP

#include <cstdint>
#include <string>

enum Nucleotide { A = 0, C = 1, G = 2, T = 3, N = 4 };

constexpr Nucleotide char_to_nucleotide(char c);

class Kmer {
  public:
    using data_t = std::uint64_t;
    explicit Kmer(std::size_t K) : _data(0), K(K), n_count(0) {}
    explicit Kmer(const std::string &kmer);
    void roll(char c);

    /**
     * @brief Get the nucleotide at position i (0-indexed, from right to left)
     * @param i The position of the nucleotide to get
     * @return The nucleotide at position i
     */
    Nucleotide get(std::size_t i) const;
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
    data_t reverse_complement(data_t data) const;
    std::size_t K;
    data_t _data;
    std::size_t n_count;
};

#endif
