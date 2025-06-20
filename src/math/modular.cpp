#include "math/modular.hpp"
#include <algorithm>
#include <bit>
#include <cassert>

size_t bit_size(auto a) {
    size_t bits = 0;
    while (a) {
        bits++;
        a >>= 1;
    }
    return bits;
}

bool check_nooverflow(auto a, auto b, size_t bits) {
    std::size_t bitsa = bit_size(a);
    std::size_t bitsb = bit_size(b);
    return bitsa + bitsb < bits;
}

bool check_nounderflow(std::uint64_t a, std::uint64_t b) { return a >= b; }

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
    assert(check_nooverflow(x, magic, sizeof(uint128_t) * 8));
    std::uint64_t quotient = (((uint128_t)magic * x) >> exp);
    assert(check_nooverflow(quotient, mod, sizeof(uint128_t) * 8));
    assert(check_nounderflow(x, quotient * mod));
    assert(bit_size(x - (uint128_t)quotient * mod) <=
           sizeof(std::uint64_t) * 8);
    std::uint64_t remainder = x - (uint128_t)quotient * mod;
    return remainder;
}
