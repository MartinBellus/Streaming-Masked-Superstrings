#include "bloom_filter/rolling_hash.hpp"

std::uint64_t rolling_hash::operator()(const std::string &key) const {
    std::uint64_t hash = 0;
    for (char c : key) {
        hash = mod.reduce2(hash * p + c);
    }
    return mod.reduce(hash);
}

constexpr std::uint64_t primes[] = {31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
constexpr std::uint64_t mods[] = {
        1000000000000099669, 1000000000000099729, 1000000000000099733,
        1000000000000099739, 1000000000000099759, 1000000000000099853,
        1000000000000099873, 1000000000000099907, 1000000000000099909,
        1000000000000099937, 1000000000000099957, 1000000000000099961};

rolling_hash_family::rolling_hash_family(std::size_t nhashes)
    : nhashes(nhashes) {
    hash_fn.reserve(nhashes);
    for (std::size_t i = 0; i < nhashes; i++) {
        hash_fn.emplace_back(primes[i], mods[i]);
    }
}

std::vector<std::uint64_t>
rolling_hash_family::operator()(const std::string &key) const {
    std::vector<std::uint64_t> hashes;
    hashes.reserve(nhashes);
    for (auto &&h : hash_fn) {
        hashes.push_back(h(key));
    }
    return hashes;
}
