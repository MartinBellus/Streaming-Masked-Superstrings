#ifndef HASHED_KMER_HPP
#define HASHED_KMER_HPP

#include "hash/hash_family.hpp"
#include "helper/kmer.hpp"

template <class T>
class HashedKmer;

template <Hash H>
class HashedKmer<H> : public Kmer {
  public:
    HashedKmer(std::size_t K, std::uint64_t seed, KmerRepr repr)
        : Kmer(K), _hash(seed, repr) {}
    H::hash_t hash() const { return _hash.hash(*this); }

  private:
    H _hash;
};

template <RollingHash H>
class HashedKmer<H> : public Kmer {
  public:
    HashedKmer(std::size_t K, std::uint64_t seed, KmerRepr repr)
        : Kmer(K), repr(repr), _hash(K, seed) {}
    void roll(char c) {
        Nucleotide n_in = char_to_nucleotide(c);
        Nucleotide n_out = Kmer::last(KmerRepr::FORWARD);
        _hash.roll(n_in, n_out);

        Kmer::roll(c);
    }
    void reset() {
        Kmer::reset();
        _hash.reset();
    }
    H::hash_t hash() const {
        bool use_reverse = Kmer::use_reverse(repr);
        return _hash.get_hash(use_reverse);
    }

  private:
    KmerRepr repr;
    H _hash;
};

struct KmerHash {
    template <class H>
    std::size_t operator()(const HashedKmer<H> &kmer) const {
        return kmer.hash();
    }
};

#endif
