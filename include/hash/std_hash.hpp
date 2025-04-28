#ifndef STD_HASH_HPP
#define STD_HASH_HPP

#include "hash_family.hpp"

class std_hash_family : public hash_family {
  public:
    std_hash_family(std::size_t nhashes) : hash_family(nhashes) {}
    std::span<const hash_t> hash_impl(const Kmer &key) const {
        auto h = std::hash<Kmer::data_t>();
        auto x = h(key.data());
        auto y = h(key.data() << 1);

        for (std::size_t i = 0; i < nhashes; i++) {
            buffer[i] = x + i * y;
        }
        return std::span(buffer.get(), nhashes);
    }
};

#endif
