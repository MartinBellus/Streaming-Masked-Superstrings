#ifndef COUNTING_BLOOM_FILTER_HPP
#define COUNTING_BLOOM_FILTER_HPP

#include "hash/hash_family.hpp"
#include "helper/counting_bitset.hpp"
#include "math/modular.hpp"
#include <cmath>

template <HashFamily H, std::size_t BPC = 4>
class CountingBloomFilter {
    using Self = CountingBloomFilter;

  public:
    static Self optimal(std::size_t num_elements, std::size_t bits_per_element,
                        KmerRepr repr) {
        std::size_t size = num_elements * bits_per_element;
        std::size_t nhashes = std::round(std::log(2) * bits_per_element);
        return Self(size, nhashes, repr);
    }
    CountingBloomFilter(std::size_t size, std::size_t nhashes, KmerRepr repr)
        : _size(size), hash_family(nhashes, repr), data(size) {}
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
            if (!data.is_stuck(_size.reduce(h))) {
                data.decrement(_size.reduce(h));
            }
        }
    }
    bool contains(const Kmer &key) const {
        bool contains = true;
        for (auto &&h : hash_family.hash(key)) {
            contains &= data.test(_size.reduce(h));
        }
        return contains;
    }
    std::size_t size() const { return _size.get_mod() * BPC; }
    double error_rate(std::size_t num_elements) const {
        std::size_t k = hash_family.size();
        double p = std::pow(
                1 - std::exp(-(double)k * num_elements / _size.get_mod()), k);
        return p;
    }

  private:
    CountingBitset<BPC> data;
    Modulus _size;
    mutable H hash_family;
};

template <RollingHashFamily H, std::size_t BPC = 4>
class RollingCountingBloomFilter {
    using Self = RollingCountingBloomFilter;

  public:
    static Self optimal(std::size_t num_elements, std::size_t bits_per_element,
                        std::size_t k, KmerRepr repr) {
        std::size_t size = num_elements * bits_per_element;
        std::size_t nhashes = std::round(std::log(2) * bits_per_element);
        return Self(size, nhashes, k, repr);
    }
    RollingCountingBloomFilter(std::size_t size, std::size_t nhashes,
                               std::size_t k, KmerRepr repr)
        : _size(size), hash_family(nhashes, k, repr), data(size) {}
    void init(const Kmer &key) { hash_family.init(key); }
    void reset_hash_family() { hash_family.reset(); }
    void roll(char c) { hash_family.roll(c); }
    void insert_this() {
        if (contains_this()) {
            return;
        }
        for (auto &&h : hash_family.get_hashes()) {
            data.increment(_size.reduce(h));
        }
    }
    void erase_this() {
        if (!contains_this()) {
            return;
        }
        for (auto &&h : hash_family.get_hashes()) {
            if (!data.is_stuck(_size.reduce(h))) {
                data.decrement(_size.reduce(h));
            }
        }
    }
    bool contains_this() const {
        bool contains = true;
        for (auto &&h : hash_family.get_hashes()) {
            contains &= data.test(_size.reduce(h));
        }
        return contains;
    }
    static std::size_t bucket_size() { return BPC; }
    std::size_t size() const { return _size.get_mod(); }
    double error_rate(std::size_t num_elements) const {
        std::size_t k = hash_family.size();
        double p = std::pow(
                1 - std::exp(-(double)k * num_elements / _size.get_mod()), k);
        return p;
    }

  private:
    CountingBitset<BPC> data;
    Modulus _size;
    H hash_family;
};

#endif
