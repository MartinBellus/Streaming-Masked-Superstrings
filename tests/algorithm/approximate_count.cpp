#include "algorithm/approximate_count.hpp"
#include "hash/murmur_hash.hpp"
#include "hash/poly_hash.hpp"
#include "helper/args.hpp"
#include "helper/hashed_kmer.hpp"
#include "io/fasta.hpp"
#include <chrono>
#include <iostream>
#include <unordered_set>

std::size_t exact_count(io::FastaReader &in, std::size_t K) {
    using kmer_t = HashedKmer<poly_hash>;
    in.reset();
    std::unordered_set<kmer_t, KmerHash> kmer_set;
    kmer_t kmer(K, 0, KmerRepr::FORWARD);
    char c;
    while (in.next_sequence()) {
        kmer.reset();
        for (std::size_t i = 0; i < K && in.next_nucleotide(c); i++) {
            kmer.roll(c);
        }
        while (in.next_nucleotide(c)) {
            kmer.roll(c);
            kmer_set.insert(kmer);
        }
    }
    return kmer_set.size();
}

std::size_t exact_duplicated_count(io::FastaReader &in, std::size_t K) {
    using kmer_t = HashedKmer<poly_hash>;
    in.reset();
    std::unordered_set<kmer_t, KmerHash> kmer_set;
    std::unordered_set<kmer_t, KmerHash> duplicated_set;
    kmer_t kmer(K, 0, KmerRepr::FORWARD);
    char c;
    while (in.next_sequence()) {
        kmer.reset();
        for (std::size_t i = 0; i < K && in.next_nucleotide(c); i++) {
            kmer.roll(c);
        }
        while (in.next_nucleotide(c)) {
            kmer.roll(c);
            if (!kmer_set.insert(kmer).second) {
                duplicated_set.insert(kmer);
            }
        }
    }
    return duplicated_set.size();
}

template <typename H>
void stats_test(const std::string &path, std::size_t K) {
    io::FastaReader in(path);
    std::string fake_args[] = {"-k", std::to_string(K), path, ""};
    auto arg = ComputeArgs::from_cmdline(4, fake_args).value();
    auto start = std::chrono::high_resolution_clock::now();
    auto stats = approximate_count<H>(arg);
    auto end = std::chrono::high_resolution_clock::now();

    std::int64_t count_real = exact_count(in, K);
    std::int64_t count_duplicated = exact_duplicated_count(in, K);
    std::int64_t approximate_duplicated = stats.total_length -
                                          stats.sequence_count * (K - 1) -
                                          stats.approximate_kmer_count;
    double count_error =
            static_cast<double>(count_real - stats.approximate_kmer_count) /
            static_cast<double>(count_real) * 100;
    double duplicated_error =
            static_cast<double>(
                    std::abs(count_duplicated - approximate_duplicated)) /
            static_cast<double>(count_duplicated) * 100;

    std::cout << "Time: \t\t\t"
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                         .count()
              << "ms\n";
    std::cout << "Kmer count:\n\tApproximate:\t" << stats.approximate_kmer_count
              << "\n\tReal:\t\t" << count_real << "\n\tError:\t\t"
              << count_error << "%\n";
    std::cout << "Kmer duplicate "
              << "count:\n\tApproximate:\t" << approximate_duplicated
              << "\n\tReal:\t\t" << count_duplicated << "\n\tError:\t\t"
              << duplicated_error << "%\n";
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file>\n";
        return 1;
    }

    std::string path = argv[1];
    std::size_t K = 31;
    stats_test<murmur_hash_family>(path, K);
}
