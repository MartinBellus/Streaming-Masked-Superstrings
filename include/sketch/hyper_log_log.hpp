#ifndef HYPER_LOG_LOG_HPP
#define HYPER_LOG_LOG_HPP

#include "hash/hash_family.hpp"
#include "helper/counting_bitset.hpp"
#include "helper/kmer.hpp"
#include <cmath>

template <HashFamily H, std::size_t B = 12, std::size_t BPB = 5>
class HyperLogLog {
  private:
    using Bitset = CountingBitset<BPB>;
    static constexpr std::size_t buckets = (1 << B);
    static constexpr std::size_t counter_max = (1 << BPB) - 1;
    static constexpr std::size_t index_mask = (1 << B) - 1;
    static constexpr double alpha = 0.7213 / (1 + 1.079 / buckets);

  public:
    HyperLogLog(KmerRepr repr) : hash_family(1, repr), bitset(buckets) {}
    void update(const Kmer &k) {
        auto hash = hash_family.hash(k)[0];
        auto index = hash & index_mask;
        auto zeros = std::countl_zero(hash) + 1;

        if (bitset.get(index) < zeros) {
            bitset.set(index, zeros);
        }
    }
    std::size_t query() const {
        auto hll_est = hll_estimate();
        if (2 * hll_est <= 5 * buckets) {
            return linear_estimate();
        } else {
            return hll_est;
        }
    }

  private:
    std::size_t hll_estimate() const {
        std::size_t sum = 0;
        for (std::size_t i = 0; i < buckets; i++) {
            sum += 1 << (counter_max - bitset.get(i));
        }
        std::size_t base = static_cast<std::size_t>(alpha * buckets * buckets);
        return (base << counter_max) / sum;
    }
    std::size_t linear_estimate() const {
        std::size_t zero_count = 0;
        for (std::size_t i = 0; i < buckets; i++) {
            if (bitset.get(i) == 0) {
                zero_count++;
            }
        }
        if (zero_count == 0) {
            return buckets * std::log(buckets);
        }
        return static_cast<std::size_t>(
                buckets * std::log(static_cast<double>(buckets) / zero_count));
    }
    H hash_family;
    Bitset bitset;
};

#endif
