#ifndef BLOOM_FILTER_HPP
#define BLOOM_FILTER_HPP

#include "hash/hash_family.hpp"
#include "helper/bitset.hpp"
#include "math/modular.hpp"
#include <cmath>

template <HashFamily H>
class BloomFilter {
  public:
    static BloomFilter<H> optimal(std::size_t num_elements,
                                  std::size_t bits_per_element, KmerRepr repr) {
        std::size_t size = num_elements * bits_per_element;
        std::size_t nhashes = std::round(std::log(2) * bits_per_element);
        return BloomFilter<H>(size, nhashes, repr);
    }
    BloomFilter(std::size_t size, std::size_t nhashes, KmerRepr repr)
        : _size(size), hash_family(nhashes, repr), data(size) {}
    void insert(const Kmer &key) {
        for (auto &&h : hash_family.hash(key)) {
            data.set(_size.reduce(h));
        }
    }
    bool contains(const Kmer &key) const {
        bool contains = true;
        for (auto &&h : hash_family.hash(key)) {
            contains &= data.test(_size.reduce(h));
        }
        return contains;
    }
    std::size_t size() const { return _size.get_mod(); }
    double error_rate(std::size_t num_elements) const {
        std::size_t k = hash_family.size();
        double p = std::pow(
                1 - std::exp(-(double)k * num_elements / _size.get_mod()), k);
        return p;
    }

  private:
    DynamicBitset data;
    Modulus _size;
    mutable H hash_family;
};

template <RollingHashFamily H>
class RollingBloomFilter {
  public:
    static RollingBloomFilter<H> optimal(std::size_t num_elements,
                                         std::size_t bits_per_element,
                                         std::size_t k, KmerRepr repr) {
        std::size_t size = num_elements * bits_per_element;
        std::size_t nhashes = std::round(std::log(2) * bits_per_element);
        return RollingBloomFilter<H>(size, nhashes, k, repr);
    }
    RollingBloomFilter(std::size_t size, std::size_t nhashes, std::size_t k,
                       KmerRepr repr)
        : _size(size), hash_family(nhashes, k, repr), data(size) {}
    void init(const Kmer &key) { hash_family.init(key); }
    void reset_hash_family() { hash_family.reset(); }
    void roll(char c) { hash_family.roll(c); }
    void insert_this() {
        for (auto &&h : hash_family.get_hashes()) {
            data.set(_size.reduce(h));
        }
    }
    bool contains_this() const {
        bool contains = true;
        for (auto &&h : hash_family.get_hashes()) {
            contains &= data.test(_size.reduce(h));
        }
        return contains;
    }
    bool contains(const Kmer &kmer) const {
        H tmp_hash_family(hash_family);
        bool contains = true;
        for (auto &&h : tmp_hash_family.hash(kmer)) {
            contains &= data.test(_size.reduce(h));
        }
        return contains;
    }
    std::size_t size() const { return _size.get_mod(); }
    double error_rate(std::size_t num_elements) const {
        std::size_t k = hash_family.size();
        double p = std::pow(
                1 - std::exp(-(double)k * num_elements / _size.get_mod()), k);
        return p;
    }

  private:
    DynamicBitset data;
    Modulus _size;
    H hash_family;
};

#endif
