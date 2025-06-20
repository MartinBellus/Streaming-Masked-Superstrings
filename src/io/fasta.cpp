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
    int to_print = std::min(kmer.available(), kmer.size() - 1);
    for (int i = to_print - 1; i >= 0; i--) {
        stream.write(nucleotide_to_char[kmer.get(i)]);
    }
    kmer.reset();
}

void KmerWriter::write_header(const std::string &header) {
    stream.write('>');
    stream.write(header);
    stream.write('\n');
}
