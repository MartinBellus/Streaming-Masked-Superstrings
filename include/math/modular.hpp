#ifndef MODULAR_HPP
#define MODULAR_HPP

#include <cstdint>

using uint128_t = __uint128_t;

class Modulus {
  public:
    Modulus(std::uint64_t mod);
    std::uint64_t reduce(uint128_t x) const;
    std::uint64_t reduce2(uint128_t x) const;
    std::uint64_t get_mod() const { return mod; }

  private:
    std::uint64_t mod, magic;
};

#endif
