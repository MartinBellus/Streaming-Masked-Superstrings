#include "math/modular.hpp"
#include <algorithm>
#include <bit>

Modulus::Modulus(std::uint64_t mod) : mod(mod) {
    std::uint8_t log2 = sizeof(std::uint64_t) * 8 - std::countl_zero(mod);
    exp = std::min(60 + log2, 80);
    magic = ((uint128_t)1 << exp) / mod;
}

std::uint64_t Modulus::reduce(uint128_t x) const {
    std::uint64_t remainder = reduce2(x);
    return remainder - mod * (remainder >= mod);
}

std::uint64_t Modulus::reduce2(uint128_t x) const {
    std::uint64_t quotient = (((uint128_t)magic * x) >> exp);
    std::uint64_t remainder = x - quotient * mod;
    return remainder;
}
