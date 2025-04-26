#include "io/fasta.hpp"
#include <cctype>

using namespace io;

constexpr char nucleotide_to_char[] = {'a', 'c', 'g', 't', 'X'};

void KmerWriter::print_nucleotide(int present) {
    auto n = kmer.last();
    auto to_print = nucleotide_to_char[n];
    if (present == PRESENT) {
        to_print = std::toupper(to_print);
    }
    stream.write(to_print);
}

void KmerWriter::flush() {
    for (int i = kmer.size() - 1; i >= 0; i--) {
        stream.write(nucleotide_to_char[kmer.get(i)]);
    }
    kmer.reset();
}
