#include "helper/kmer.hpp"
#include <cassert>
#include <iostream>
#include <string>

std::uint64_t reverse(std::uint64_t data, std::size_t K) {
    data = (data >> 2 & 0x3333333333333333ULL) |
           (data << 2 & 0xCCCCCCCCCCCCCCCCULL);
    data = (data >> 4 & 0x0F0F0F0F0F0F0F0FULL) |
           (data << 4 & 0xF0F0F0F0F0F0F0F0ULL);
    data = (data >> 8 & 0x00FF00FF00FF00FFULL) |
           (data << 8 & 0xFF00FF00FF00FF00ULL);
    data = (data >> 16 & 0x0000FFFF0000FFFFULL) |
           (data << 16 & 0xFFFF0000FFFF0000ULL);
    data = (data >> 32) | (data << 32);

    auto mask = (1ULL << (2 * K)) - 1;
    data >>= 8 * sizeof(data) - 2 * K;

    return ~data & mask;
}

char nucleotide_to_char[] = {'A', 'C', 'G', 'T', 'N'};

std::string kmer_to_string(const Kmer &kmer, KmerRepr repr) {
    std::string result;
    for (int i = kmer.size() - 1; i >= 0; i--) {
        Nucleotide nucleotide = kmer.get(i, repr);
        result += nucleotide_to_char[nucleotide];
    }
    if (repr == KmerRepr::REVERSE) {
        std::reverse(result.begin(), result.end());
    }
    return result;
}

std::string reverse(const std::string &s) {
    std::string reversed;
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        char c = *it;
        switch (c) {
        case 'A':
            reversed += 'T';
            break;
        case 'C':
            reversed += 'G';
            break;
        case 'G':
            reversed += 'C';
            break;
        case 'T':
            reversed += 'A';
            break;
        }
    }
    return reversed;
}

std::string random_dna(size_t N) {
    std::string result;
    result.reserve(N);
    for (size_t i = 0; i < N; i++) {
        result += nucleotide_to_char[rand() % 4];
    }
    return result;
}

bool reverse_test(size_t N, size_t K) {
    std::string s = random_dna(N);
    Kmer kmer(K);
    size_t i = 0;
    for (; i < K - 1; i++) {
        kmer.roll(s[i]);
    }
    for (size_t start = 0; i < N; i++, start++) {
        kmer.roll(s[i]);
        auto str = kmer_to_string(kmer, KmerRepr::FORWARD);
        if (kmer.data(KmerRepr::REVERSE) !=
            reverse(kmer.data(KmerRepr::FORWARD), kmer.size())) {
            std::cerr << "Data mismatch at position " << start << ": "
                      << kmer.data(KmerRepr::REVERSE) << " != "
                      << reverse(kmer.data(KmerRepr::FORWARD), kmer.size())
                      << std::endl;
            return false;
        }
        if (str != s.substr(start, K)) {
            std::cerr << "Forward mismatch at position " << start << ": " << str
                      << " != " << s.substr(start, K) << std::endl;
            return false;
        }
        if (kmer_to_string(kmer, KmerRepr::REVERSE) != reverse(str)) {
            std::cerr << "Reverse mismatch at position " << start << ": "
                      << kmer_to_string(kmer, KmerRepr::REVERSE)
                      << " != " << reverse(str) << std::endl;
            return false;
        }
    }
    return true;
}

bool canonical_test1(size_t N, size_t K) {
    std::string target = random_dna(K);
    std::string s = target + random_dna(N - 2 * K) + reverse(target);
    Kmer kmer1(K), kmer2(K);
    int i = 0;
    for (; i < K; i++) {
        kmer1.roll(s[i]);
        kmer2.roll(s[i]);
    }
    for (; i < N; i++) {
        kmer2.roll(s[i]);
    }
    return kmer1.data(KmerRepr::CANON) == kmer2.data(KmerRepr::CANON);
}

bool canonical_test2(size_t N, size_t K) {
    std::string target = random_dna(K);
    std::string s = target + random_dna(N - 2 * K) + target;
    Kmer kmer1(K), kmer2(K);
    int i = 0;
    for (; i < K; i++) {
        kmer1.roll(s[i]);
        kmer2.roll(s[i]);
    }
    for (; i < N; i++) {
        kmer2.roll(s[i]);
    }
    return kmer1.data(KmerRepr::CANON) == kmer2.data(KmerRepr::CANON);
}

int main() {
    for (size_t i = 1; i < 31; i++) {
        assert(reverse_test(100, i));
    }
    for (size_t i = 1; i < 31; i++) {
        assert(canonical_test1(100, i));
        assert(canonical_test2(100, i));
    }
}
