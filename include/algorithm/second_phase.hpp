#ifndef SECOND_PHASE_HPP
#define SECOND_PHASE_HPP

#include "hash/hash_family.hpp"
#include "helper/args.hpp"
#include "io/fasta.hpp"
#include "sketch/counting_bloom_filter.hpp"
#include <cctype>

namespace second_phase {

template <RollingHashFamily H>
int compute_superstring(std::size_t approx_set_size, io::FastaReader &in,
                        io::KmerWriter &out, const ComputeArgs &arg) {
    using CBF = RollingCountingBloomFilter<H>;
    auto K = arg.k();
    auto repr = arg.unidirectional() ? KmerRepr::FORWARD : KmerRepr::CANON;
    CBF filter = CBF::optimal(approx_set_size, arg.bits_per_element(), repr);

    in.reset();
    while (in.next_sequence()) {
        char c;
        std::uint64_t mask = 0;
        while (in.next_nucleotide(c)) {
            filter.roll(c);
            mask = (mask << 1) | (bool)std::islower(c);
            mask &= (1ULL << K) - 1;
            if (mask & (1ULL << (K - 1))) {
                filter.insert_this();
            }
        }
    }

    in.reset();
    while (in.next_sequence()) {
        char c;
        std::uint64_t mask = 0;
        while (in.next_nucleotide(c)) {
            filter.roll(c);
            mask = (mask << 1) | (bool)std::isupper(c);
            mask &= (1ULL << K) - 1;
            if (mask & (1ULL << (K - 1))) {
                filter.erase_this();
            }
        }
    }

    in.reset();
    out.write_header(arg.fasta_header() + " (second phase)");
    while (in.next_sequence()) {
        char c;
        std::uint64_t mask = 0;
        while (in.next_nucleotide(c)) {
            filter.roll(c);
            mask = (mask << 1) | (bool)std::isupper(c);
            mask &= (1ULL << K) - 1;
            out.add_nucleotide(c);
            if (mask & (1ULL << (K - 1)) || filter.contains_this()) {
                out.print_nucleotide(io::PRESENT);
            } else {
                out.print_nucleotide(io::NOT_PRESENT);
            }
            filter.erase_this();
        }
        out.flush();
    }
}

} // namespace second_phase
#endif
