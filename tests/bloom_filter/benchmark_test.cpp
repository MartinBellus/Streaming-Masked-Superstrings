#include "bloom_filter/bloom_filter.hpp"
#include "bloom_filter/rolling_hash.hpp"
#include "bloom_filter/std_hash.hpp"
#include <cassert>
#include <chrono>
#include <iostream>
#include <unordered_set>

using namespace std;

string rand_str(size_t len, char from = 'a', char to = 'z') {
    string s(len, ' ');
    for (char &c : s) {
        c = rand() % (to - from + 1) + from;
    }
    return s;
}

template <class T>
void benchmark(size_t num, const string &name) {
    cout << endl << "Benchmarking " << name << endl;
    T bf(7);
    unordered_set<string> strings;
    for (size_t i = 0; i < num; i++) {
        strings.insert(rand_str(100));
    }
    auto insert_start = chrono::high_resolution_clock::now();
    for (auto s : strings) {
        bf.insert(s);
    }
    auto insert_end = chrono::high_resolution_clock::now();

    size_t false_positives = 0;
    vector<string> queries(num);
    vector<bool> anss(num, false);
    for (size_t i = 0; i < num; i++) {
        queries[i] = rand_str(100);
        anss[i] = strings.contains(queries[i]);
    }

    auto query_start = chrono::high_resolution_clock::now();
    for (auto &&s : strings) {
        assert(bf.contains(s));
    }
    for (size_t i = 0; i < num; i++) {
        if (bf.contains(queries[i]) != anss[i]) {
            false_positives++;
        }
    }
    auto query_end = chrono::high_resolution_clock::now();

    auto percent = (double)false_positives / (double)num * 100;
    cout << "False positive rate: " << false_positives << " / " << num << ", "
         << percent << "%" << endl;
    cout << "Insert time per element: "
         << chrono::duration_cast<chrono::nanoseconds>(insert_end -
                                                       insert_start)
                            .count() /
                    (double)num
         << " ns" << endl;
    cout << "Query time per element: "
         << chrono::duration_cast<chrono::nanoseconds>(query_end - query_start)
                            .count() /
                    (double)(2 * num)
         << " ns" << endl;
}

int main() {
    const size_t NUM = 1000000;
    using bf1 = BloomFilter<10 * NUM, rolling_hash_family>;
    using bf2 = BloomFilter<10 * NUM, std_hash_family>;

    benchmark<bf1>(NUM, "Bloom filter, rolling hash");
    benchmark<bf2>(NUM, "Bloom filter, std hash");
    benchmark<unordered_set<string>>(NUM, "Unordered set");
}
