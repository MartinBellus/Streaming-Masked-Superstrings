#include "algorithm/approximate_count.hpp"
#include "hash/murmur_hash.hpp"
#include "hash/poly_hash.hpp"
#include "helper/hashed_kmer.hpp"
#include <iostream>
#include <unordered_set>

std::size_t exact_count(io::FastaReader &in, std::size_t K) {
    using kmer_t = HashedKmer<poly_hash>;
    std::unordered_set<kmer_t, KmerHash> kmer_set;
    kmer_t kmer(K, 0);
    char c;
    while (in.next_sequence()) {
        kmer.reset();
        for (std::size_t i = 0; i < K && in.next_nucleotide(c); i++) {
            kmer.roll(c);
        }
        while (in.next_nucleotide(c)) {
            auto nout = kmer.last();
            kmer.roll(c);
            kmer_set.insert(kmer);
        }
    }
    return kmer_set.size();
}

template <typename H>
void hll_test(const std::string &path, std::size_t K) {
    io::FastaReader in(path);
    std::int64_t count = approximate_count<H>(in, K);
    in.reset();
    std::int64_t count_real = exact_count(in, K);

    std::cout << "Approximate count: \t" << count << "\n";
    std::cout << "Exact count: \t\t" << count_real << "\n";
    std::cout << "Error: "
              << static_cast<double>(count_real - count) /
                         static_cast<double>(count_real) * 100
              << "%\n";
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file>\n";
        return 1;
    }

    std::string path = argv[1];
    std::size_t K = 31;
    hll_test<murmur_hash_family>(path, K);
}
