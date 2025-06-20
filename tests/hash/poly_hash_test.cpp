#include "hash/poly_hash.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

constexpr char nucleotide_to_char[] = {'A', 'C', 'G', 'T'};

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

bool test_poly_hash_init(size_t K) {
    poly_hash hash1(K, 0), hash2(K, 0);
    std::string s = random_dna(K);
    std::string rev_s = reverse(s);

    hash1.init(Kmer(s));
    hash2.init(Kmer(rev_s));

    bool res = hash1.get_hash(false) == hash2.get_hash(true) &&
               hash1.get_hash(true) == hash2.get_hash(false);
    if (!res) {
        std::cerr << "Init test failed for K = " << K << ", s = " << s
                  << ", rev_s = " << rev_s << "\n";
        std::cerr << "hash1: " << hash1.get_hash(false) << ", "
                  << hash1.get_hash(true) << "\n";
        std::cerr << "hash2: " << hash2.get_hash(false) << ", "
                  << hash2.get_hash(true) << "\n";
    }
    return res;
}

bool test_poly_hash_roll_simple(size_t K) {
    poly_hash hash1(K, 0), hash2(K, 0);
    std::string s = random_dna(K);
    std::string rev_s = reverse(s);

    for (int i = 0; i < K; i++) {
        Nucleotide n = char_to_nucleotide(s[i]);
        Nucleotide rev_n = char_to_nucleotide(rev_s[i]);
        hash1.roll(n, N);
        hash2.roll(rev_n, N);
    }

    bool res = hash1.get_hash(false) == hash2.get_hash(true) &&
               hash1.get_hash(true) == hash2.get_hash(false);
    if (!res) {
        std::cerr << "Roll simple test failed for K = " << K << ", s = " << s
                  << ", rev_s = " << rev_s << "\n";
        std::cerr << "hash1: " << hash1.get_hash(false) << ", "
                  << hash1.get_hash(true) << "\n";
        std::cerr << "hash2: " << hash2.get_hash(false) << ", "
                  << hash2.get_hash(true) << "\n";
    }
    return res;
}

bool test_poly_hash_roll(size_t n, size_t K) {
    std::string target = random_dna(K);
    std::string s = random_dna(n - K) + target;
    std::string rev_s = random_dna(n - K) + reverse(target);
    poly_hash hash1(K, 0), hash2(K, 0);
    int i = 0;
    for (; i < K; i++) {
        auto n_in = char_to_nucleotide(s[i]);
        auto rn_in = char_to_nucleotide(rev_s[i]);
        hash1.roll(n_in, N);
        hash2.roll(rn_in, N);
    }
    for (; i < n; i++) {
        auto n_in = char_to_nucleotide(s[i]);
        auto n_out = char_to_nucleotide(s[i - K]);
        auto rn_in = char_to_nucleotide(rev_s[i]);
        auto rn_out = char_to_nucleotide(rev_s[i - K]);
        hash1.roll(n_in, n_out);
        hash2.roll(rn_in, rn_out);
    }
    bool res = hash1.get_hash(false) == hash2.get_hash(true) &&
               hash1.get_hash(true) == hash2.get_hash(false);
    if (!res) {
        std::cerr << "Roll test failed for K = " << K << ", s = " << s
                  << ", rev_s = " << rev_s << "\n";
        std::cerr << "hash1: " << hash1.get_hash(false) << ", "
                  << hash1.get_hash(true) << "\n";
        std::cerr << "hash2: " << hash2.get_hash(false) << ", "
                  << hash2.get_hash(true) << "\n";
    }
    return res;
}

bool test_poly_hash_complements(size_t n, size_t K) {
    std::string target = random_dna(K);
    std::string s = random_dna(n) + target + random_dna(n) + reverse(target);
    poly_hash hash(K, 0);
    int i = 0;
    for (; i < K; i++) {
        auto n_in = char_to_nucleotide(s[i]);
        hash.roll(n_in, N);
    }
    for (; i < n + K; i++) {
        auto n_in = char_to_nucleotide(s[i]);
        auto n_out = char_to_nucleotide(s[i - K]);
        hash.roll(n_in, n_out);
    }
    auto res1 = hash.get_hash(false);
    for (; i < s.size(); i++) {
        auto n_in = char_to_nucleotide(s[i]);
        auto n_out = char_to_nucleotide(s[i - K]);
        hash.roll(n_in, n_out);
    }
    auto res2 = hash.get_hash(true);
    bool res = res1 == res2;
    if (!res) {
        std::cerr << "Complement test failed for K = " << K << ", n = " << n
                  << ", s = " << s << "\n";
        std::cerr << "hash: " << res1 << ", " << res2 << "\n";
    }
    return res;
}

bool test_hash_family(size_t N, size_t K) {
    std::string target = random_dna(K);
    std::string s =
            target + random_dna(std::min((size_t)0, N - K)) + reverse(target);
    poly_hash_family hashF(4, K, KmerRepr::CANON);

    size_t i = 0;
    for (; i < K; i++) {
        hashF.roll(s[i]);
    }

    auto &&hash1 = hashF.get_hashes();
    std::vector<std::uint64_t> hashes1(hash1.begin(), hash1.end());

    for (; i < s.size(); i++) {
        hashF.roll(s[i]);
    }

    auto &&hash2 = hashF.get_hashes();
    std::vector<std::uint64_t> hashes2(hash2.begin(), hash2.end());

    bool res = hashes1 == hashes2;
    if (!res) {
        std::cerr << "Hash family test failed for K = " << K << ", N = " << N
                  << ", s = " << s << "\n";
        std::cerr << "hash1: ";
        for (const auto &h : hashes1) {
            std::cerr << h << " ";
        }
        std::cerr << "\nhash2: ";
        for (const auto &h : hashes2) {
            std::cerr << h << " ";
        }
        std::cerr << "\n";
    }
    return res;
}

int main() {
    for (int k = 1; k < 31; k++) {
        for (int i = 0; i < 100; i++) {
            assert(test_poly_hash_init(k));
        }
    }

    for (int k = 1; k < 31; k++) {
        for (int i = 0; i < 100; i++) {
            assert(test_poly_hash_roll_simple(k));
        }
    }

    for (int k = 1; k < 31; k++) {
        for (int i = 0; i < 100; i++) {
            assert(test_poly_hash_roll(100, k));
        }
    }

    for (int k = 1; k < 31; k++) {
        for (int i = 0; i < 50; i++) {
            assert(test_poly_hash_complements(k / 2, k));
        }
        for (int i = 0; i < 50; i++) {
            assert(test_poly_hash_complements(2 * k, k));
        }
    }

    for (int k = 1; k < 31; k++) {
        for (int i = 0; i < 50; i++) {
            assert(test_hash_family(k / 2, k));
        }
        for (int i = 0; i < 50; i++) {
            assert(test_hash_family(2 * k, k));
        }
    }
}
