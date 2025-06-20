#ifndef KMER_HPP
#define KMER_HPP

#include <cstdint>
#include <string>

enum Nucleotide { A = 0, C = 1, G = 2, T = 3, N = 4 };
constexpr Nucleotide COMPLEMENT[] = {T, G, C, A, N};

Nucleotide char_to_nucleotide(char c);

enum KmerRepr {
    FORWARD,
    REVERSE,
    CANON,
};

class Kmer {
  public:
    using data_t = std::uint64_t;
    explicit Kmer(std::size_t K) : K(K), _data(0), _rev_data(0), n_count(0) {}
    explicit Kmer(const std::string &kmer);
    void roll(char c);

    /**
     * @brief Get the nucleotide at position i (0-indexed, from right to left)
     * @param i The position of the nucleotide to get
     * @param representation Whether to use forward, reverse or canonical kmer
     * representation
     * @return The nucleotide at position i
     */
    Nucleotide get(std::size_t i, KmerRepr representation) const;
    Nucleotide last(KmerRepr representation) const {
        return get(K - 1, representation);
    }
    std::size_t size() const { return K; }
    std::size_t available() const { return std::min(K, n_count); }
    void reset() {
        _data = 0;
        _rev_data = 0;
        n_count = 0;
    }
    bool operator==(const Kmer &other) const {
        return other.K == K &&
               (_data == other._data || _data == other._rev_data);
    }
    const data_t &data(KmerRepr representation) const;
    bool use_reverse(KmerRepr representation) const;

  private:
    std::size_t K;
    data_t _data, _rev_data;
    std::size_t n_count;
};

#endif
