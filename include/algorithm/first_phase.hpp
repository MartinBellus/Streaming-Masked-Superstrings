#ifndef FIRST_PHASE_HPP
#define FIRST_PHASE_HPP

#include "helper/args.hpp"
#include "io/fasta.hpp"

namespace first_phase {

template <class BF>
int compute_superstring(std::size_t approx_set_size, io::FastaReader &in,
                        io::KmerWriter &out, const ComputeArgs &args) {
    auto K = args.k();
    in.reset();
    out.write_header(args.fasta_header());
    BF filter(10 * approx_set_size, 7, K);
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
