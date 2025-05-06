#include "hash/murmur_hash.hpp"

constexpr std::uint64_t m = 0xc6a4a7935bd1e995ULL;
constexpr std::uint32_t r = 47;
constexpr std::uint64_t mask = 0xff;

constexpr char nucleotide_to_char[] = {'a', 'c', 'g', 't', 'X'};

std::string kmer_to_string(const Kmer &kmer) {
    std::string s;
    s.reserve(kmer.size());
    for (std::size_t i = 0; i < kmer.size(); i++) {
        s.push_back(nucleotide_to_char[kmer.get(i)]);
    }
    return s;
}

murmur_hash::hash_t murmur_hash::hash(const Kmer &key) const {
    auto key_str = kmer_to_string(key);
    auto len = key_str.size();
    hash_t h = _seed ^ (len * m);

    const uint64_t *data = (const uint64_t *)key_str.data();
    const uint64_t *end = data + (len / 8);

    while (data != end) {
        uint64_t k = *data++;

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const unsigned char *data2 = (const unsigned char *)data;

    switch (len & 7) {
    case 7:
        h ^= ((uint64_t)data2[6]) << 48;
    case 6:
        h ^= ((uint64_t)data2[5]) << 40;
    case 5:
        h ^= ((uint64_t)data2[4]) << 32;
    case 4:
        h ^= ((uint64_t)data2[3]) << 24;
    case 3:
        h ^= ((uint64_t)data2[2]) << 16;
    case 2:
        h ^= ((uint64_t)data2[1]) << 8;
    case 1:
        h ^= ((uint64_t)data2[0]);
        h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

murmur_hash_family::murmur_hash_family(std::size_t nhashes)
    : hash_family(nhashes), xhash(0x1234), yhash(0x5678) {}

std::span<const murmur_hash_family::hash_t>
murmur_hash_family::hash_impl(const Kmer &kmer) {
    auto x = xhash.hash(kmer);
    auto y = yhash.hash(kmer);
    for (std::size_t i = 0; i < nhashes; i++) {
        buffer[i] = x + i * y;
    }
    return std::span(buffer.get(), nhashes);
}
