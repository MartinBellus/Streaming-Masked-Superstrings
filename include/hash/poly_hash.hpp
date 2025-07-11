#ifndef ROLLING_HASH_HPP
#define ROLLING_HASH_HPP

#include "hash_family.hpp"
#include "helper/kmer.hpp"
#include "math/modular.hpp"
#include <cstdint>

class poly_hash {
  public:
    static constexpr bool rolling = true;
    using hash_t = std::uint64_t;
    poly_hash(std::size_t k, std::uint64_t seed);
    poly_hash(std::uint64_t p, std::uint64_t mod, std::size_t k);
    poly_hash(std::uint64_t p, std::uint64_t mod, const Kmer &kmer)
        : p(p), mod(mod) {
        init(kmer);
    }
    hash_t get_hash(bool reverse) const { return reverse ? rev_state : state; };
    void roll(Nucleotide n_in, Nucleotide n_out);
    void init(const Kmer &kmer);
    void reset();

  private:
    std::uint64_t p, inv_p, state, rev_state, last_exp;
    std::size_t k;
    Modulus mod;
};

class poly_hash_family : public rolling_hash_family<poly_hash_family> {
  public:
    poly_hash_family(std::size_t nhashes, std::size_t k, KmerRepr repr);
    poly_hash_family(std::size_t nhashes, KmerRepr repr);
    void roll_impl(char c);
    void init_impl(const Kmer &kmer);
    void reset_impl();

  private:
    void update_hashes();
    poly_hash xhash, yhash;
    KmerRepr repr;
    Kmer kmer;
};

static_assert(RollingHash<poly_hash>);
static_assert(RollingHashFamily<poly_hash_family>);

#endif
