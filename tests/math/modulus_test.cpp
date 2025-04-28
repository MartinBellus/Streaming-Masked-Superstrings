#include "math/modular.hpp"
#include <iostream>
#include <random>
#include <stdexcept>

using namespace std;

mt19937_64 rng;

void test_all(uint64_t mod, uint64_t lo, uint64_t hi, uint64_t step = 1) {
    Modulus m(mod);
    for (uint64_t i = lo; i < hi; i += step) {
        auto x = m.reduce(i);
        if (x != i % mod) {
            throw runtime_error("Modulus test failed: " + to_string(i) + " % " +
                                to_string(mod) + " != " + to_string(x) +
                                ", expected " + to_string(i % mod));
        }
    }
}

__uint128_t random_128(__uint128_t lo, __uint128_t hi) {
    auto rnd = ((__uint128_t)rng() << 64) + rng();
    return rnd % (hi - lo) + lo;
}

void test_random(uint64_t mod, __uint128_t lo, __uint128_t hi, size_t count) {
    Modulus m(mod);
    for (size_t i = 0; i < count; i++) {
        auto x = random_128(lo, hi);
        auto a = m.reduce(x);
        if (a != x % mod) {
            auto res = x % mod;
            uint64_t res_high = res >> 64;
            uint64_t res_low = res;
            uint64_t high = x >> 64;
            uint64_t low = x;
            throw runtime_error("Modulus test failed: " + to_string(high) +
                                to_string(low) + " % " + to_string(mod) +
                                " != " + to_string(a) + ", expected " +
                                to_string(res_high) + to_string(res_low));
        }
    }
}

int main() {
    uniform_int_distribution<uint64_t> dist;
    dist = uniform_int_distribution<uint64_t>(1, 10000);
    for (size_t i = 0; i < 100; i++) {
        auto mod = dist(rng);
        test_all(mod, 0, 10 * mod);
    }
    cerr << "Small modulus with small numbers OK" << endl;

    for (size_t m = 1; m < 40; m++) {
        dist = uniform_int_distribution<uint64_t>(1ULL << m, 1ULL << (m + 1));
        for (size_t i = 0; i < 64; i++) {
            auto mod = dist(rng);
            test_random(mod, (__uint128_t)1 << i, (__uint128_t)1 << (i + 1),
                        1ULL << (m / 2 + 1));
        }
    }
    cerr << "All bit combinations OK" << endl;

    dist = uniform_int_distribution<uint64_t>(1e5, 1e9);
    for (size_t i = 0; i < 100; i++) {
        auto mod = dist(rng);
        test_all(mod, 0, 10 * mod, 1033);
    }
    cerr << "Medium modulus with medium numbers OK" << endl;

    dist = uniform_int_distribution<uint64_t>(1ULL << 40, 1ULL << 50);
    for (size_t i = 0; i < 10; i++) {
        auto mod = dist(rng);
        test_random(mod, mod, 10 * mod, (int)1e6);
    }
    cerr << "Big modulus with big numbers OK" << endl;
}
