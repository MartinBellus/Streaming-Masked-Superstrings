#ifndef COUNTING_BLOOM_FILTER_HPP
#define COUNTING_BLOOM_FILTER_HPP

#include "hash/hash_family.hpp"
#include "helper/counting_bitset.hpp"
#include "math/modular.hpp"

template <HashFamily H, std::size_t BPC = 4>
class CountingBloomFilter {
  public:
    CountingBloomFilter(std::size_t size, std::size_t nhashes)
        : _size(size), hash_family(nhashes), data(size) {}
    void insert(const Kmer &key) {
        if (contains(key)) {
            return;
        }
        for (auto &&h : hash_family.hash(key)) {
            data.increment(_size.reduce(h));
        }
    }
    void erase(const Kmer &key) {
        if (!contains(key)) {
            return;
        }
        for (auto &&h : hash_family.hash(key)) {
            data.decrement(_size.reduce(h));
        }
    }
    bool contains(const Kmer &key) const {
        bool contains = true;
        for (auto &&h : hash_family.hash(key)) {
            contains &= data.test(_size.reduce(h));
        }
        return contains;
    }

  private:
    CountingBitset<BPC> data;
    Modulus _size;
    mutable H hash_family;
};

template <RollingHashFamily H, std::size_t BPC = 4>
class RollingCountingBloomFilter {
  public:
    RollingCountingBloomFilter(std::size_t size, std::size_t nhashes,
                               std::size_t k)
        : _size(size), hash_family(nhashes, k), data(size) {}
    void init(const Kmer &key) { hash_family.init(key); }
    void reset_hash_family() { hash_family.reset(); }
    void roll(char c) { hash_family.roll(c); }
    void insert_this() {
        if (contains_this()) {
            return;
        }
        for (auto &&h : hash_family.get_hashes()) {
            data.inc(_size.reduce(h));
        }
    }
    void erase_this() {
        if (!contains_this()) {
            return;
        }
        for (auto &&h : hash_family.get_hashes()) {
            data.dec(_size.reduce(h));
        }
    }
    bool contains_this() const {
        bool contains = true;
        for (auto &&h : hash_family.get_hashes()) {
            contains &= data.test(_size.reduce(h));
        }
        return contains;
    }

  private:
    CountingBitset<BPC> data;
    Modulus _size;
    H hash_family;
};

#endif
