#include "algorithm/exact.hpp"
#include "hash/poly_hash.hpp"
#include "helper/hashed_kmer.hpp"
#include <unordered_set>

using namespace exact;

int exact::compute_superstring(std::size_t K, io::FastaReader &in,
                               io::KmerWriter &out) {
    using kmer_t = HashedKmer<poly_hash>;
    std::unordered_set<kmer_t, KmerHash> kmer_set;
    kmer_t kmer(K, 0);
    char c;
    while (in.next_sequence()) {
        kmer.reset();
        for (std::size_t i = 0; i < K - 1 && in.next_nucleotide(c); i++) {
            kmer.roll(c);
            out.add_nucleotide(c);
        }
        while (in.next_nucleotide(c)) {
            kmer.roll(c);
            out.add_nucleotide(c);
            if (kmer_set.contains(kmer)) {
                out.print_nucleotide(io::NOT_PRESENT);
            } else {
                out.print_nucleotide(io::PRESENT);
                kmer_set.insert(kmer);
            }
        }
        out.flush();
    }
    return 0;
}
