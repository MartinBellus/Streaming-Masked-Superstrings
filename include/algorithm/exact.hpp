#ifndef EXACT_HPP
#define EXACT_HPP

#include "helper/args.hpp"

namespace exact {

struct Accuracy {
    std::size_t missing_kmers = 0;
    std::size_t additional_kmers = 0;
    std::size_t length = 0;
};

Accuracy compute_accuracy(const CompareArgs &args);

int compute_superstring(const ExactArgs &args);
} // namespace exact

std::ostream &operator<<(std::ostream &os, const exact::Accuracy &acc);

#endif
