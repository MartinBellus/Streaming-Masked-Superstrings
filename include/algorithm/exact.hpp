#ifndef EXACT_HPP
#define EXACT_HPP

#include "io/fasta.hpp"

namespace exact {
int compute_superstring(std::size_t K, io::FastaReader &in,
                        io::KmerWriter &out);
}

#endif
