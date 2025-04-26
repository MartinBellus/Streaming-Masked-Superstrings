#ifndef ROLLING_HASH_HPP
#define ROLLING_HASH_HPP

#include "hash_family.hpp"
#include "helper/kmer.hpp"
#include "math/modular.hpp"
#include <cstdint>
#include <string>

class poly_hash {
  public:
    poly_hash(std::uint64_t p, std::uint64_t mod, std::size_t k);
    poly_hash(std::uint64_t p, std::uint64_t mod, const Kmer &kmer)
        : p(p), mod(mod) {
        init(kmer);
    }
    std::uint64_t get_hash() const { return state; };
    void roll(Nucleotide n_in, Nucleotide n_out);
    void init(const Kmer &kmer);
    void reset();

  private:
    std::uint64_t p, state, last_exp;
    std::size_t k;
    Modulus mod;
};

class poly_hash_family : public rolling_hash_family {
  public:
    poly_hash_family(std::size_t nhashes, std::size_t k);
    poly_hash_family(std::size_t nhashes);
    void roll_impl(char c);
    void init_impl(const std::string &s);
    void reset_impl();

  private:
    poly_hash xhash, yhash;
    Kmer kmer;
};

#endif
