#include "hash/poly_hash.hpp"
#include "helper/kmer.hpp"

constexpr std::uint64_t NVALUE[] = {1, 2, 3, 4};

constexpr std::uint64_t primes[] = {31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
// constexpr std::uint64_t mods[] = {
//         1000000000000099669, 1000000000000099729, 1000000000000099733,
//         1000000000000099739, 1000000000000099759, 1000000000000099853,
//         1000000000000099873, 1000000000000099907, 1000000000000099909,
//         1000000000000099937, 1000000000000099957, 1000000000000099961};
constexpr std::uint64_t mods[] = {1099511627781, 1099511627789, 1099511627793,
                                  1099511627799, 1099511627807, 1099511627811,
                                  1099511627819, 1099511627823, 1099511627829,
                                  1099511627843};

std::uint64_t pow_mod(std::uint64_t a, std::uint64_t b, const Modulus &mod) {
    std::uint64_t result = 1;
    while (b) {
        if (b % 2 == 1) {
            result = mod.reduce2((uint128_t)result * a);
        }
        a = mod.reduce2((uint128_t)a * a);
        b /= 2;
    }
    return mod.reduce(result);
}

poly_hash::poly_hash(std::size_t k)
    : p(primes[0]), mod(mods[0]), state(0), k(k) {
    last_exp = pow_mod(p, k, mod);
}

poly_hash::poly_hash(std::uint64_t p, std::uint64_t mod, std::size_t k)
    : p(p), mod(mod), k(k) {
    last_exp = pow_mod(p, k, mod);
}

void poly_hash::init(const Kmer &kmer) {
    state = 0;
    k = kmer.size();
    last_exp = pow_mod(p, k, mod);
    for (std::size_t i = 0; i < kmer.size(); i++) {
        auto nucleotide = kmer.get(kmer.size() - i - 1);
        state = mod.reduce2((uint128_t)state * p + NVALUE[nucleotide]);
    }
    state = mod.reduce(state);
}
void poly_hash::roll(Nucleotide n_in, Nucleotide n_out) {
    state = mod.reduce2((uint128_t)state * p + NVALUE[n_in]);
    auto last = mod.reduce2(NVALUE[n_out] * last_exp);
    state = mod.reduce2(2 * mod.get_mod() + state - last);
    state = mod.reduce(state);
}

void poly_hash::reset() { state = 0; }

poly_hash_family::poly_hash_family(std::size_t nhashes, std::size_t k)
    : rolling_hash_family(nhashes), kmer(k), xhash(primes[0], mods[0], k),
      yhash(primes[1], mods[1], k) {}

poly_hash_family::poly_hash_family(std::size_t nhashes)
    : poly_hash_family(nhashes, 0) {}

void poly_hash_family::roll_impl(char c) {
    Nucleotide n_in = char_to_nucleotide(c);
    Nucleotide n_out = kmer.get(kmer.size() - 1);
    kmer.roll(c);
    xhash.roll(n_in, n_out);
    yhash.roll(n_in, n_out);
    auto x = xhash.get_hash();
    auto y = yhash.get_hash();
    for (std::size_t i = 0; i < nhashes; i++) {
        buffer[i] = x + i * y;
    }
}

void poly_hash_family::init_impl(const Kmer &key) {
    kmer = key;
    xhash.init(kmer);
    yhash.init(kmer);
    auto x = xhash.get_hash();
    auto y = yhash.get_hash();
    for (std::size_t i = 0; i < nhashes; i++) {
        buffer[i] = x + i * y;
    }
}

void poly_hash_family::reset_impl() {
    xhash.reset();
    yhash.reset();
    kmer.reset();
}
