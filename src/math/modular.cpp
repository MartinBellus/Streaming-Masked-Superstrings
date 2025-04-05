#include "math/modular.hpp"

Modulus::Modulus(std::uint64_t mod) : mod(mod), magic((-1ULL) / mod) {}

std::uint64_t Modulus::reduce(uint128_t x) const {
    std::uint64_t remainder = reduce2(x);
    return remainder - mod * (remainder >= mod);
}

std::uint64_t Modulus::reduce2(uint128_t x) const {
    std::uint64_t quotient = (((uint128_t)magic * x) >> 64);
    std::uint64_t remainder = x - quotient * mod;
    return remainder;
}
