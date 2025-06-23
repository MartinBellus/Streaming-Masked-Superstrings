#ifndef APPROXIMATE_COUNT_HPP
#define APPROXIMATE_COUNT_HPP

#include "hash/hash_family.hpp"
#include "helper/args.hpp"
#include "io/fasta.hpp"
#include "sketch/hyper_log_log.hpp"

template <HashFamily H>
std::size_t approximate_count(const ComputeArgs &arg) {
    io::FastaReader in(arg.dataset());
    auto K = arg.k();
    auto kmer_repr = arg.unidirectional() ? KmerRepr::FORWARD : KmerRepr::CANON;
    HyperLogLog<H> hll(kmer_repr);
    Kmer kmer(K);
    char c;

    while (in.next_sequence()) {
        kmer.reset();
        for (std::size_t i = 0; i < K && in.next_nucleotide(c); i++) {
            kmer.roll(c);
        }
        while (in.next_nucleotide(c)) {
            hll.update(kmer);
            kmer.roll(c);
        }
    }

    return hll.query();
}

#endif
