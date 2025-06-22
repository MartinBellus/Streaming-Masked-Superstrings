#include "io/fasta.hpp"
#include <cctype>

using namespace io;

constexpr char nucleotide_to_char[] = {'a', 'c', 'g', 't', 'X'};

void KmerWriter::print_nucleotide(int present) {
    auto n = kmer.last();
    auto to_print = nucleotide_to_char[n];
    if (present == PRESENT) {
        last_one = 0;
        to_print = std::toupper(to_print);
    }
    if (!splice || last_one < kmer.size()) {
        stream.write(to_print);
    }
    last_one++;
}

void KmerWriter::flush() {
    int to_print = std::min(kmer.available(), kmer.size() - 1);
    for (int i = to_print - 1; i >= 0; i--) {
        add_nucleotide('A');
        print_nucleotide(NOT_PRESENT);
    }
    kmer.reset();
}

void KmerWriter::write_header(const std::string &header) {
    stream.write('>');
    stream.write(header);
    stream.write('\n');
}
