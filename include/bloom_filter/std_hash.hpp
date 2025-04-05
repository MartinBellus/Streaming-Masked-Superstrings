#ifndef STD_HASH_HPP
#define STD_HASH_HPP

#include "hash_family.hpp"
#include <string>

class std_hash_family : public hash_family {
  public:
    std_hash_family(std::size_t nhashes) : hash_family(nhashes) {}
    std::span<const hash_t> hash_impl(const std::string &key) const {
        auto h = std::hash<std::string>();
        auto x = h(key);
        auto y = h(key + '0');

        for (std::size_t i = 0; i < nhashes; i++) {
            buffer[i] = x + i * y;
        }
        return std::span(buffer.get(), nhashes);
    }
};

#endif
