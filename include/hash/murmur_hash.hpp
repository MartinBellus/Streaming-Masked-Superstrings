#ifndef MURMUR_HASH_HPP
#define MURMUR_HASH_HPP

#include "hash_family.hpp"

class murmur_hash {
  public:
    static constexpr bool rolling = false;
    using hash_t = std::uint64_t;
    murmur_hash(std::uint64_t seed) : _seed(seed) {}
    hash_t hash(const Kmer &key) const;
    std::uint64_t seed() const { return _seed; }

  private:
    std::uint64_t _seed;
};

class murmur_hash_family : public hash_family<murmur_hash_family> {

  public:
    murmur_hash_family(std::size_t nhashes);
    std::span<const hash_t> hash_impl(const Kmer &kmer);

  private:
    murmur_hash xhash, yhash;
};

static_assert(Hash<murmur_hash>);
static_assert(HashFamily<murmur_hash_family>);

#endif
