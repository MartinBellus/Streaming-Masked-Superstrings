#include "helper/fasta.hpp"
#include <fstream>
#include <iostream>

int readfile(const std::string &name) {
    std::ifstream file(name);
    if (!file) {
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
