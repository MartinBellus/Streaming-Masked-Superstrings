#ifndef ROLLING_HASH_HPP
#define ROLLING_HASH_HPP

#include "hash_family.hpp"
#include "math/modular.hpp"
#include <cstdint>
#include <string>
#include <vector>

class rolling_hash {
  public:
    rolling_hash(std::uint64_t p, std::uint64_t mod) : p(p), mod(mod) {}
    std::uint64_t operator()(const std::string &key) const;

  private:
    std::uint64_t p;
    Modulus mod;
};

class rolling_hash_family : public hash_family_tag {
  public:
    rolling_hash_family(std::size_t nhashes);
    std::vector<std::uint64_t> operator()(const std::string &key) const;

  private:
    std::size_t nhashes;
    std::vector<rolling_hash> hash_fn;
};

#endif
