#ifndef STD_HASH_HPP
#define STD_HASH_HPP

#include "hash_family.hpp"
#include <cstdint>
#include <string>
#include <vector>

class std_hash_family : public hash_family_tag {
  public:
    std_hash_family(std::size_t nhashes) : nhashes(nhashes) {}
    std::vector<std::uint64_t> operator()(const std::string &key) const {
        std::vector<std::uint64_t> hashes;
        hashes.reserve(nhashes);
        auto h = std::hash<std::string>();
        auto x = h(key);
        auto y = h(key + '0');

        for (std::size_t i = 0; i < nhashes; i++) {
            hashes.push_back((x + i * y));
        }
        return hashes;
    }

  private:
    std::size_t nhashes;
};

#endif
