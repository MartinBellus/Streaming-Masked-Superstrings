#include "algorithm/approximate_count.hpp"
#include "algorithm/exact.hpp"
#include "algorithm/first_phase.hpp"
#include "algorithm/second_phase.hpp"
#include "hash/murmur_hash.hpp"
#include "hash/poly_hash.hpp"
#include "helper/args.hpp"
#include <iostream>
#include <ranges>
#include <vector>

using namespace io;

int usage() {
    // clang-format off
    std::cerr << std::endl;
    std::cerr << "Program: StreamingMaskedSuperstrings - a tool for computation of approximate masked superstrings with small memory footprint." << std::endl;
    std::cerr << "Usage:   streaming-masked-superstrings <command> [options]" << std::endl << std::endl;
    std::cerr << "Command:" << std::endl;
    std::cerr << "    compute    - Compute an approximate masked superstring from input FASTA." << std::endl;
    std::cerr << "    exact      - Compute an exact masked superstring from input FASTA." << std::endl;
    std::cerr << "    compare    - Compare exact and approximate masked superstrings and report accuracy." << std::endl;
    std::cerr << std::endl;
    // clang-format on
    return 1;
}

int subcomand_compute(auto &&args) {
    auto _arg = ComputeArgs::from_cmdline(args.size(), args.data());
    if (!_arg.has_value()) {
        return ComputeArgs::usage();
    }
    auto arg = _arg.value();
    auto stats = approximate_count<murmur_hash_family>(arg);

    auto ret = first_phase::compute_superstring<poly_hash_family>(
            stats.approximate_kmer_count, arg);

    if (arg.second_phase()) {
        auto approximate_duplicates = stats.total_length -
                                      stats.sequence_count * (arg.k() - 1) -
                                      stats.approximate_kmer_count;
        return second_phase::compute_superstring<poly_hash_family>(
                approximate_duplicates, arg);
    }

    return ret;
}

int subcomand_exact(auto &&args) {
    auto _arg = ExactArgs::from_cmdline(args.size(), args.data());
    if (!_arg.has_value()) {
        return ExactArgs::usage();
    }
    auto arg = _arg.value();
    return exact::compute_superstring(arg);
}

int subcomand_compare(auto &&args) {
    auto _arg = CompareArgs::from_cmdline(args.size(), args.data());
    if (!_arg.has_value()) {
        return CompareArgs::usage();
    }
    auto arg = _arg.value();
    auto acc = exact::compute_accuracy(arg);
    std::cout << acc;
    return 0;
}

int main(int argc, char *argv[]) {
    std::vector<std::string> args(argv, argv + argc);
    if (argc == 1 || args[1] == "--help" || args[1] == "-h") {
        return usage();
    }
    auto subcommand_args = args | std::views::drop(2);
    if (args[1] == "compute") {
        return subcomand_compute(subcommand_args);
    }
    if (args[1] == "exact") {
        return subcomand_exact(subcommand_args);
    }
    if (args[1] == "compare") {
        return subcomand_compare(subcommand_args);
    }
    return usage();
}
