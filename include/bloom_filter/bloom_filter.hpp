#ifndef BLOOM_FILTER_HPP
#define BLOOM_FILTER_HPP

#include "hash_family.hpp"
#include <bitset>
#include <string>

template <std::size_t S, HashFamily H>
class BloomFilter {
  public:
    BloomFilter(std::size_t nhashes) : hash_family(nhashes) {}
    void insert(const std::string &key) {
        for (auto &&h : hash_family.hash(key)) {
            data.set((h % S));
        }
    }
    bool contains(const std::string &key) const {
        bool contains = true;
        for (auto &&h : hash_family.hash(key)) {
            contains &= data.test(h % S);
        }
        return contains;
    }

  private:
    std::bitset<S> data;
    mutable H hash_family;
};

template <std::size_t S, RollingHashFamily H>
class RollingBloomFilter {
  public:
    RollingBloomFilter(std::size_t nhashes, std::size_t k)
        : hash_family(nhashes, k) {}
    void init(const std::string &key) { hash_family.init(key); }
    void roll(char c) { hash_family.roll(c); }
    void insert_this() {
        for (auto &&h : hash_family.get_hashes()) {
            data.set(h % S);
        }
    }
    bool contains_this() const {
        bool contains = true;
        for (auto &&h : hash_family.get_hashes()) {
            contains &= data.test(h % S);
        }
        return contains;
    }
    bool contains(const std::string &s) {
        bool contains = true;
        for (auto &&h : hash_family.hash(s)) {
            contains &= data.test(h % S);
        }
        return contains;
    }

  private:
    std::bitset<S> data;
    H hash_family;
};

#endif
