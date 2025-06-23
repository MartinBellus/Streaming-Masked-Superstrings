#ifndef FIRST_PHASE_HPP
#define FIRST_PHASE_HPP

#include "hash/hash_family.hpp"
#include "helper/args.hpp"
#include "io/fasta.hpp"
#include "sketch/bloom_filter.hpp"

namespace first_phase {

template <RollingHashFamily H>
int compute_superstring(std::size_t approx_set_size, const ComputeArgs &args) {
    using BF = RollingBloomFilter<H>;
    auto K = args.k();
    bool splice = args.splice() && !args.second_phase();
    io::FastaReader in(args.dataset());
    io::KmerWriter out(args.first_phase_output(), K, splice);
    auto kmer_repr =
            args.unidirectional() ? KmerRepr::FORWARD : KmerRepr::CANON;
    out.write_header(args.fasta_header());
    BF filter =
            BF::optimal(approx_set_size, args.bits_per_element(), K, kmer_repr);
    while (in.next_sequence()) {
        std::size_t read = 0;
        char c;
        filter.reset_hash_family();
        for (std::size_t i = 0; i < K - 1 && in.next_nucleotide(c); i++) {
            filter.roll(c);
            out.add_nucleotide(c);
        }
        while (in.next_nucleotide(c)) {
            filter.roll(c);
            out.add_nucleotide(c);
            bool first_occurence = !filter.contains_this();
            if (first_occurence) {
                filter.insert_this();
                out.print_nucleotide(io::PRESENT);
            } else {
                out.print_nucleotide(io::NOT_PRESENT);
            }
        }
        out.flush();
    }
    return 0;
}
} // namespace first_phase

#endif
