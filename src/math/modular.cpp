#include "math/modular.hpp"

Modulus::Modulus(std::uint64_t mod) : mod(mod), magic((1ULL << 32) / mod) {}

using u128 = __uint128_t;

std::uint64_t Modulus::reduce(std::uint64_t x) const {
    std::uint64_t quotient = (((u128)magic * x) >> 64);
    std::uint64_t remainder = x - quotient * magic;
    return remainder - mod * (remainder >= mod);
}

std::uint64_t Modulus::reduce2(std::uint64_t x) const {
    std::uint64_t quotient = (((u128)magic * x) >> 64);
    std::uint64_t remainder = x - quotient * magic;
    return remainder;
}
