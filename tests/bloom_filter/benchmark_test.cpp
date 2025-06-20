#include "hash/murmur_hash.hpp"
#include "hash/poly_hash.hpp"
#include "sketch/bloom_filter.hpp"
#include <chrono>
#include <iostream>
#include <unordered_set>
#include <vector>

using namespace std;

constexpr string ALPHABET = "ACGT";

string rand_seq(size_t len) {
    string s(len, ' ');
    for (char &c : s) {
        c = ALPHABET[rand() % ALPHABET.size()];
    }
    return s;
}

template <class T>
void roll_benchmark(size_t len, size_t k, const string &name) {
    cout << endl << "Benchmarking " << name << endl;
    T bf(10 * len, 7, k, KmerRepr::FORWARD);
    string s = rand_seq(len + k - 1);
    unordered_set<string> strings;
    for (size_t i = 0; i < len; i++) {
        strings.insert(s.substr(i, k));
    }
    auto insert_start = chrono::high_resolution_clock::now();
    bf.init(Kmer(s.substr(0, k)));
    bf.insert_this();
    for (size_t i = k; i < s.size(); i++) {
        bf.roll(s[i]);
        bf.insert_this();
    }
    auto insert_end = chrono::high_resolution_clock::now();

    size_t false_positives = 0;
    vector<string> queries(len);
    vector<bool> anss(len, false);
    for (size_t i = 0; i < len; i++) {
        queries[i] = rand_seq(k);
        anss[i] = strings.contains(queries[i]);
    }

    auto query_start = chrono::high_resolution_clock::now();
    for (auto &&s : strings) {
        if (!bf.contains(Kmer(s))) {
            throw runtime_error("Filter " + name + " does not contain " + s);
        }
    }
    for (size_t i = 0; i < len; i++) {
        if (bf.contains(Kmer(queries[i])) != anss[i]) {
            false_positives++;
        }
    }
    auto query_end = chrono::high_resolution_clock::now();

    auto percent = (double)false_positives / (double)len * 100;
    cout << "False positive rate: " << false_positives << " / " << len << ", "
         << percent << "%" << endl;
    cout << "Insert time per element: "
         << chrono::duration_cast<chrono::nanoseconds>(insert_end -
                                                       insert_start)
                            .count() /
                    (double)len
         << " ns" << endl;
    cout << "Query time per element (non sequential): "
         << chrono::duration_cast<chrono::nanoseconds>(query_end - query_start)
                            .count() /
                    (double)(2 * len)
         << " ns" << endl;
}

template <typename T>
concept is_unordered_set = std::is_same_v<T, std::unordered_set<std::string>>;

template <class T>
T create(std::size_t size, std::size_t nhashes) {
    if constexpr (is_unordered_set<T>) {
        return T();
    } else {
        return T(size, nhashes, KmerRepr::FORWARD);
    }
}

template <class T>
void benchmark(size_t len, size_t k, const string &name) {
    using key = std::conditional_t<is_unordered_set<T>, std::string, Kmer>;

    cout << endl << "Benchmarking " << name << endl;
    T bf = create<T>(10 * len, 7);
    unordered_set<string> strings;
    for (size_t i = 0; i < len; i++) {
        strings.insert(rand_seq(k));
    }
    auto insert_start = chrono::high_resolution_clock::now();
    for (auto &&s : strings) {
        bf.insert(key(s));
    }
    auto insert_end = chrono::high_resolution_clock::now();

    size_t false_positives = 0;
    vector<string> queries(len);
    vector<bool> anss(len, false);
    for (size_t i = 0; i < len; i++) {
        queries[i] = rand_seq(k);
        anss[i] = strings.contains(queries[i]);
    }

    auto query_start = chrono::high_resolution_clock::now();
    for (auto &&s : strings) {
        if (!bf.contains(key(s))) {
            throw runtime_error("Filter " + name + " does not contain " + s);
        }
    }
    for (size_t i = 0; i < len; i++) {
        bool res;
        if (bf.contains(key(queries[i])) != anss[i]) {
            false_positives++;
        }
    }
    auto query_end = chrono::high_resolution_clock::now();

    auto percent = (double)false_positives / (double)len * 100;
    cout << "False positive rate: " << false_positives << " / " << len << ", "
         << percent << "%" << endl;
    cout << "Insert time per element: "
         << chrono::duration_cast<chrono::nanoseconds>(insert_end -
                                                       insert_start)
                            .count() /
                    (double)len
         << " ns" << endl;
    cout << "Query time per element: "
         << chrono::duration_cast<chrono::nanoseconds>(query_end - query_start)
                            .count() /
                    (double)(2 * len)
         << " ns" << endl;
}

int main() {
    const size_t NUM = 1000000;
    const size_t K = 31;
    using bf1 = BloomFilter<poly_hash_family>;
    using bf2 = BloomFilter<murmur_hash_family>;
    using rbf1 = RollingBloomFilter<poly_hash_family>;

    benchmark<bf1>(NUM, K, "Bloom filter, rolling hash");
    benchmark<bf2>(NUM, K, "Bloom filter, murmur hash");
    benchmark<unordered_set<string>>(NUM, K, "Unordered set");
    roll_benchmark<rbf1>(NUM, K, "Rolling bloom filter, rolling hash");
}
