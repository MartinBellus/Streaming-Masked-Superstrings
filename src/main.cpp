#include "algorithm/approximate_count.hpp"
#include "algorithm/exact.hpp"
#include "algorithm/first_phase.hpp"
#include "hash/murmur_hash.hpp"
#include "hash/poly_hash.hpp"
#include "io/fasta.hpp"
#include "sketch/bloom_filter.hpp"
#include <iostream>
#include <vector>

using namespace io;
using rolling_bf = RollingBloomFilter<poly_hash_family>;
constexpr std::size_t K = 31;

int usage() {
    std::cerr << "Usage: superstring [--exact] <fasta-input> <output>"
              << std::endl;
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        return usage();
    }
    std::vector<std::string> args(argv, argv + argc);
    if (argc == 4) {
        if (args[1] == "--exact") {
            FastaReader in(argv[2]);
            KmerWriter out(argv[3], K);
            return exact::compute_superstring(K, in, out);
        } else if (args[1] == "--compare") {
            FastaReader output(argv[2]);
            FastaReader golden_output(argv[3]);
            auto acc = exact::compute_accuracy(output, golden_output);
            std::cout << acc;
            return 0;
        } else {
            return usage();
        }
    }
    FastaReader in(argv[1]);
    KmerWriter out(argv[2], K);
    auto size = approximate_count<murmur_hash_family>(in, K);
    return first_phase::compute_superstring<rolling_bf>(K, size, in, out);
}
