#ifndef HYPER_LOG_LOG_HPP
#define HYPER_LOG_LOG_HPP

#include "hash/hash_family.hpp"
#include "helper/counting_bitset.hpp"
#include "helper/kmer.hpp"

template <HashFamily H, std::size_t BPB = 5, std::size_t B = 3>
class HyperLogLog {
  private:
    using Bitset = CountingBitset<BPB>;
    using inner_t = typename Bitset::inner_t;
    static constexpr std::size_t inner_size = sizeof(inner_t) * 8;
    static constexpr std::size_t index_mask = (1 << B) - 1;

  public:
    HyperLogLog() : hash_family(1), bitset(1 << B) {}
    void update(const Kmer &k) {
        auto hash = hash_family.hash(k)[0];
        auto index = hash & index_mask;
        auto zeros = std::countl_zero(hash);

        if (bitset.get(index) < zeros) {
            bitset.set(index, zeros);
        }
    }
    std::size_t query() const {
        std::size_t sum = 0;
        for (std::size_t i = 0; i < bitset.size(); i++) {
            sum += 1 << (inner_size - bitset.get(i));
        }
        return (1ULL << inner_size) / sum;
    }

  private:
    H hash_family;
    Bitset bitset;
};

#endif
