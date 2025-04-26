#ifndef BLOOM_FILTER_HPP
#define BLOOM_FILTER_HPP

#include "hash_family.hpp"
#include "helper/bitset.hpp"
#include "math/modular.hpp"
#include <string>

template <HashFamily H>
class BloomFilter {
  public:
    BloomFilter(std::size_t size, std::size_t nhashes)
        : _size(size), hash_family(nhashes), data(size) {}
    void insert(const std::string &key) {
        for (auto &&h : hash_family.hash(key)) {
            data.set(_size.reduce(h));
        }
    }
    bool contains(const std::string &key) const {
        bool contains = true;
        for (auto &&h : hash_family.hash(key)) {
            contains &= data.test(_size.reduce(h));
        }
        return contains;
    }

  private:
    DynamicBitset data;
    Modulus _size;
    mutable H hash_family;
};

template <RollingHashFamily H>
class RollingBloomFilter {
  public:
    RollingBloomFilter(std::size_t size, std::size_t nhashes, std::size_t k)
        : _size(size), hash_family(nhashes, k), data(size) {}
    void init(const std::string &key) { hash_family.init(key); }
    void reset_hash_family() { hash_family.reset(); }
    void roll(char c) { hash_family.roll(c); }
    void insert_this() {
        for (auto &&h : hash_family.get_hashes()) {
            data.set(_size.reduce2(h));
        }
    }
    bool contains_this() const {
        bool contains = true;
        for (auto &&h : hash_family.get_hashes()) {
            contains &= data.test(_size.reduce2(h));
        }
        return contains;
    }
    bool contains(const std::string &s) {
        bool contains = true;
        for (auto &&h : hash_family.hash(s)) {
            contains &= data.test(_size.reduce2(h));
        }
        return contains;
    }

  private:
    DynamicBitset data;
    Modulus _size;
    H hash_family;
};

#endif
