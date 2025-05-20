#ifndef EXACT_HPP
#define EXACT_HPP

#include "helper/args.hpp"
#include "io/fasta.hpp"

namespace exact {

struct Accuracy {
    std::size_t missing_kmers = 0;
    std::size_t additional_kmers = 0;
    std::size_t length = 0;
};

Accuracy compute_accuracy(io::FastaReader &output,
                          io::FastaReader &golden_output);

int compute_superstring(io::FastaReader &in, io::KmerWriter &out,
                        const ExactArgs &args);
} // namespace exact

std::ostream &operator<<(std::ostream &os, const exact::Accuracy &acc);

#endif
