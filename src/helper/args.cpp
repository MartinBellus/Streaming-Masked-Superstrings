#include "helper/args.hpp"
#include <filesystem>
#include <format>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

using opt_set = std::unordered_set<std::string>;

bool next_opt(std::string *&begin, const std::string *end, std::string &opt,
              std::string &val, const opt_set &opts, const opt_set &flags) {
    if (begin == end) {
        return false;
    }
    for (auto &&f : flags) {
        if (*begin != f) {
            continue;
        }
        opt = f;
        val = "";
        begin++;
        return true;
    }
    for (auto &&o : opts) {
        if (!begin->starts_with(o)) {
            continue;
        }
        opt = o;
        if (begin->size() == o.size()) {
            if (++begin == end) {
                return false;
            }
            val = *begin;
        } else {
            val = begin->substr(o.size());
        }
        begin++;
        return true;
    }
    return false;
}

std::string get_tmp_file_name(const std::string &input) {
    std::string tmp_file_filename = std::format(
            "{}-{:05}.tmp", std::filesystem::path(input).stem().string(),
            time(NULL) % 100000);
    return std::filesystem::temp_directory_path().string() + tmp_file_filename;
}

std::optional<ComputeArgs> ComputeArgs::from_cmdline(int argc,
                                                     std::string *argv) {
    const opt_set opts = {"-k", "-bpk", "-t"};
    const opt_set flags = {"-u", "-s", "--no-splice", "-f"};
    std::unordered_map<std::string, std::string> opt_map = {{"-k", "31"},
                                                            {"-bpk", "10"}};
    auto begin = argv;
    auto end = argv + argc;
    std::string opt, val;
    while (next_opt(begin, end, opt, val, opts, flags)) {
        opt_map[opt] = val;
    }
    if (end - begin != 2) {
        return std::nullopt;
    }
    auto input = *(begin++);
    auto first_out = *(begin++);
    std::string second_out = "";
    if (!opt_map.contains("-f")) {
        if (opt_map.contains("-t")) {
            second_out = opt_map.at("-t");
        } else {
            second_out = get_tmp_file_name(input);
        }
        swap(first_out, second_out);
    }
    try {
        return ComputeArgs(
                std::stoul(opt_map.at("-k")), std::stoul(opt_map.at("-bpk")),
                opt_map.contains("-u"),
                opt_map.contains("-s") || opt_map.contains("--no-splice"),
                opt_map.contains("-f"), std::move(input), std::move(first_out),
                std::move(second_out));
    } catch (...) {
        return std::nullopt;
    }
}

int ComputeArgs::usage() {
    // clang-format off
    std::cerr << "Usage: streaming-masked-superstrings compute [options] <input-fasta> <output-fasta>" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -k <int>         kmer size (default = 31)" << std::endl;
    std::cerr << "  -bpk <int>       bits per kmer (default = 10)" << std::endl;
    std::cerr << "  -t <path>        path to the temporary file used in second phase" << std::endl;
    std::cerr << "  -u               treat kmer and its reverse complement as distinct" << std::endl;
    std::cerr << "  -s, --no-splice  do not splice the resulting masked superstring" << std::endl;
    std::cerr << "  -f               run only the first phase of the algorithm" << std::endl;
    // clang-format on
    return 1;
}

std::string ComputeArgs::fasta_header() const {
    std::string mode = _unidirectional ? "unidirectional" : "bidirectional";
    return std::format("approximate masked superstring dataset='{}' k={} "
                       "bits-per-kmer={} mode={} splice={}",
                       _dataset, _k, _bpk, mode, !_no_splice);
}

std::optional<ExactArgs> ExactArgs::from_cmdline(int argc, std::string *argv) {
    const opt_set opts = {"-k"};
    const opt_set flags = {"-u", "-s", "--no-splice"};
    std::unordered_map<std::string, std::string> opt_map = {{"-k", "31"}};
    auto begin = argv;
    auto end = argv + argc;
    std::string opt, val;
    while (next_opt(begin, end, opt, val, opts, flags)) {
        opt_map[opt] = val;
    }
    if (end - begin != 2) {
        return std::nullopt;
    }
    auto input = *(begin++);
    auto output = *(begin++);
    try {
        return ExactArgs(std::stoul(opt_map.at("-k")), opt_map.contains("-u"),
                         opt_map.contains("-s") ||
                                 opt_map.contains("--no-splice"),
                         std::move(input), std::move(output));
    } catch (...) {
        return std::nullopt;
    }
}

int ExactArgs::usage() {
    // clang-format off
    std::cerr << "Usage: streaming-masked-superstrings exact [options] <input-fasta> <output-fasta>" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -k <int>         kmer size (default = 31)" << std::endl;
    std::cerr << "  -u               treat kmer and its reverse complement as distinct" << std::endl;
    std::cerr << "  -s, --no-splice  do not splice the resulting masked superstring" << std::endl;
    // clang-format on
    return 1;
}

std::string ExactArgs::fasta_header() const {
    std::string mode = _unidirectional ? "unidirectional" : "bidirectional";
    return std::format(
            "exact masked superstring dataset='{}' k={} mode={} splice={}",
            _dataset, _k, mode, !_no_splice);
}

std::optional<CompareArgs> CompareArgs::from_cmdline(int argc,
                                                     std::string *argv) {
    if (argc != 2) {
        return std::nullopt;
    }
    return CompareArgs(std::move(argv[0]), std::move(argv[1]));
}

int CompareArgs::usage() {
    // clang-format off
    std::cerr << "Usage: streaming-masked-superstrings compare [options] <approximate-output> <exact-output>" << std::endl;
    std::cerr << "Options:" << std::endl;
    // clang-format on
    return 1;
}
