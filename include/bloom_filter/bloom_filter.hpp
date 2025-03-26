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
        for (auto &&h : hash_family(key)) {
            data.set((h % S));
        }
    }
    bool contains(const std::string &key) const {
        bool contains = true;
        for (auto &&h : hash_family(key)) {
            contains &= data.test(h % S);
        }
        return contains;
    }

  private:
    std::bitset<S> data;
    H hash_family;
};

#endif
