#ifndef APPROXIMATE_COUNT_HPP
#define APPROXIMATE_COUNT_HPP

#include "hash/hash_family.hpp"
#include "helper/args.hpp"
#include "io/fasta.hpp"
#include "sketch/hyper_log_log.hpp"

struct Stats {
    std::size_t approximate_kmer_count = 0;
    std::size_t sequence_count = 0;
    std::size_t total_length = 0;
};

template <HashFamily H>
Stats approximate_count(const ComputeArgs &arg) {
    Stats stats;
    io::FastaReader in(arg.dataset());
    auto K = arg.k();
    auto kmer_repr = arg.unidirectional() ? KmerRepr::FORWARD : KmerRepr::CANON;
    HyperLogLog<H> hll(kmer_repr);
    Kmer kmer(K);
    char c;

    while (in.next_sequence()) {
        stats.sequence_count++;
        kmer.reset();
        for (std::size_t i = 0; i < K && in.next_nucleotide(c); i++) {
            stats.total_length++;
            kmer.roll(c);
        }
        while (in.next_nucleotide(c)) {
            stats.total_length++;
            hll.update(kmer);
            kmer.roll(c);
        }
    }

    stats.approximate_kmer_count = hll.query();

    return stats;
}

#endif
