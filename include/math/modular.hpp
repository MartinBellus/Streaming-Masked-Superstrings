#ifndef MODULAR_HPP
#define MODULAR_HPP

#include <cstdint>

class Modulus {
  public:
    Modulus(std::uint64_t mod);
    std::uint64_t reduce(std::uint64_t x) const;
    std::uint64_t reduce2(std::uint64_t x) const;

  private:
    std::uint64_t mod, magic;
};

#endif
