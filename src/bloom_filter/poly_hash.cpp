#include "bloom_filter/poly_hash.hpp"
#include "helper/kmer.hpp"

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

poly_hash::poly_hash(std::uint64_t p, std::uint64_t mod, std::size_t k)
    : p(p), mod(mod), k(k) {
    last_exp = pow_mod(p, k, mod);
}

void poly_hash::init(const Kmer &kmer) {
    state = 0;
    k = kmer.size();
    last_exp = pow_mod(p, k, mod);
    for (std::size_t i = 0; i < kmer.size(); i++) {
        state = mod.reduce2((uint128_t)state * p +
                            kmer.get(kmer.size() - i - 1));
    }
    state = mod.reduce(state);
}
void poly_hash::roll(Nucleotide n_in, Nucleotide n_out) {
    state = mod.reduce2((uint128_t)state * p + n_in);
    auto last = mod.reduce2(n_out * last_exp);
    state = mod.reduce2(2 * mod.get_mod() + state - last);
    state = mod.reduce(state);
}

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
// constexpr std::uint64_t mods[] = {
// 1000000007, 1000000009, 1000000021, 1000000033, 1000000087, 1000000093,
// 1000000097, 1000000103, 1000000123, 1000000181, 1000000207, 1000000223,
// 1000000241, 1000000271, 1000000289, 1000000297};

poly_hash_family::poly_hash_family(std::size_t nhashes, std::size_t k)
    : rolling_hash_family(nhashes), kmer(k) {
    hash_fn.reserve(nhashes);
    for (std::size_t i = 0; i < nhashes; i++) {
        hash_fn.emplace_back(primes[i], mods[i], kmer.size());
    }
}

poly_hash_family::poly_hash_family(std::size_t nhashes)
    : rolling_hash_family(nhashes), kmer(0) {
    hash_fn.reserve(nhashes);
    for (std::size_t i = 0; i < nhashes; i++) {
        hash_fn.emplace_back(primes[i], mods[i], kmer.size());
    }
}

void poly_hash_family::roll_impl(char c) {
    Nucleotide n_in = char_to_nucleotide(c);
    Nucleotide n_out = kmer.get(kmer.size() - 1);
    kmer.roll(c);
    for (std::size_t i = 0; i < hash_fn.size(); i++) {
        hash_fn[i].roll(n_in, n_out);
        buffer[i] = hash_fn[i].get_hash();
    }
}

void poly_hash_family::init_impl(const std::string &s) {
    kmer = Kmer(s);
    for (std::size_t i = 0; i < hash_fn.size(); i++) {
        hash_fn[i].init(kmer);
        buffer[i] = hash_fn[i].get_hash();
    }
}
