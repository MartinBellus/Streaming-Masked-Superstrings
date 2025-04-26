#include "algorithm/first_phase.hpp"
#include "bloom_filter/bloom_filter.hpp"
#include "bloom_filter/poly_hash.hpp"
#include "io/fasta.hpp"
#include <iostream>

using namespace io;
using rolling_bf = RollingBloomFilter<poly_hash_family>;
constexpr std::size_t K = 31;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage " << argv[0] << " <fasta-input> <output>"
                  << std::endl;
        return 1;
    }
    FastaReader in(argv[1]);
    KmerWriter out(argv[2], K);
    return first_phase::compute_superstring<rolling_bf>(K, (std::size_t)2e6, in,
                                                        out);
}
