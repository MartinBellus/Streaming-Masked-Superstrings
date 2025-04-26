#include "io/fasta.hpp"
#include "io/streams.hpp"
#include <iostream>

using namespace io;

int readfile(const std::string &name) {
    io::input_stream file(name);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << name << std::endl;
        return 1;
    }

    FastaReader fasta(std::move(file));
    char nucleotide;
    while (fasta.next_sequence()) {
        std::cout << '>' << fasta.get_header() << std::endl;
        while (fasta.next_nucleotide(nucleotide)) {
            std::cout << nucleotide;
        }
        std::cout << std::endl;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    return readfile(argv[1]);
}
