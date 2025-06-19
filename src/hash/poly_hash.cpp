#include "hash/poly_hash.hpp"
#include "helper/kmer.hpp"

constexpr std::uint64_t primes[] = {31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73};
// constexpr std::uint64_t mods[] = {
//         1000000000000099669, 1000000000000099729, 1000000000000099733,
//         1000000000000099739, 1000000000000099759, 1000000000000099853,
//         1000000000000099873, 1000000000000099907, 1000000000000099909,
//         1000000000000099937, 1000000000000099957, 1000000000000099961};
constexpr std::uint64_t mods[] = {1099000000003, 1099000000037, 1099000000061,
                                  1099000000069, 1099000000183, 1099000000213,
                                  1099000000271, 1099000000289, 1099000000321,
                                  1099000000327, 1099000000411, 1099000000519,
                                  1099000000541, 1099000000543, 1099000000559,
                                  1099000000573, 1099000000577, 1099000000589};

constexpr std::size_t NUM_PRIMES = sizeof(primes) / sizeof(primes[0]);
constexpr std::size_t NUM_MODS = sizeof(mods) / sizeof(mods[0]);

std::uint64_t pow_mod(std::uint64_t a, std::uint64_t b, const Modulus &mod) {
    std::uint64_t result = 1;
    while (b) {
        if (b % 2 == 1) {
            result = mod.reduce((uint128_t)result * a);
        }
        a = mod.reduce((uint128_t)a * a);
        b /= 2;
    }
    return mod.reduce(result);
}

poly_hash::poly_hash(std::size_t k, std::uint64_t seed)
    : poly_hash(primes[seed % NUM_PRIMES], mods[seed % NUM_MODS], k) {}

poly_hash::poly_hash(std::uint64_t p, std::uint64_t _mod, std::size_t k)
    : p(p), mod(_mod), k(k) {
    last_exp = pow_mod(p, k, mod);
    inv_p = pow_mod(p, _mod - 2, mod);
    reset();
}

constexpr std::uint64_t NVALUE[] = {1, 2, 3, 4, 0};

void poly_hash::init(const Kmer &kmer) {
    reset();
    k = kmer.size();
    last_exp = pow_mod(p, k, mod);
    inv_p = pow_mod(p, mod.get_mod() - 2, mod);
    for (int i = kmer.size() - 1; i >= 0; i--) {
        auto nucleotide = kmer.get(i, KmerRepr::FORWARD);
        auto rnucleotide = kmer.get(i, KmerRepr::REVERSE);
        state = mod.reduce2((uint128_t)state * p + NVALUE[nucleotide]);
        rev_state = mod.reduce2((uint128_t)rev_state * p + NVALUE[rnucleotide]);
    }
    state = mod.reduce(state);
    rev_state = mod.reduce(rev_state);
}

void poly_hash::roll(Nucleotide n_in, Nucleotide n_out) {
    state = mod.reduce2((uint128_t)state * p + NVALUE[n_in]);
    auto last = mod.reduce2(NVALUE[n_out] * last_exp);
    state = mod.reduce2(2 * mod.get_mod() + state - last);
    state = mod.reduce(state);

    n_in = COMPLEMENT[n_in];
    n_out = COMPLEMENT[n_out];
    rev_state = mod.reduce2(2 * mod.get_mod() + rev_state - NVALUE[n_out] +
                            NVALUE[n_in] * last_exp);
    rev_state = mod.reduce2((uint128_t)rev_state * inv_p);
    rev_state = mod.reduce(rev_state);
}

void poly_hash::reset() {
    state = 0;
    rev_state = 0;
}

poly_hash_family::poly_hash_family(std::size_t nhashes, std::size_t k)
    : rolling_hash_family(nhashes), kmer(k), xhash(k, 0), yhash(k, 1) {}

poly_hash_family::poly_hash_family(std::size_t nhashes)
    : poly_hash_family(nhashes, 0) {}

void poly_hash_family::roll_impl(char c) {
    Nucleotide n_in = char_to_nucleotide(c);
    Nucleotide n_out = kmer.last(KmerRepr::FORWARD);
    kmer.roll(c);
    xhash.roll(n_in, n_out);
    yhash.roll(n_in, n_out);
    auto x = xhash.get_hash(0); // TODO
    auto y = yhash.get_hash(0); // TODO
    for (std::size_t i = 0; i < nhashes; i++) {
        buffer[i] = x + i * y;
    }
}

void poly_hash_family::init_impl(const Kmer &key) {
    kmer = key;
    xhash.init(kmer);
    yhash.init(kmer);
    auto x = xhash.get_hash(0); // TODO
    auto y = yhash.get_hash(0); // TODO
    for (std::size_t i = 0; i < nhashes; i++) {
        buffer[i] = x + i * y;
    }
}

void poly_hash_family::reset_impl() {
    xhash.reset();
    yhash.reset();
    kmer.reset();
}
