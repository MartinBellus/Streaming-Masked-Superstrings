#include "helper/args.hpp"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
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
    std::stringstream ss;
    ss << std::filesystem::path(input).stem().string() << "-" << std::setw(5)
       << std::setfill('0') << (time(NULL) % 100000) << ".tmp";
    std::string tmp_file_filename = ss.str();
    return std::filesystem::temp_directory_path().string() + "/" +
           tmp_file_filename;
}

std::optional<ComputeArgs> ComputeArgs::from_cmdline(int argc,
                                                     std::string *argv) {
    const opt_set opts = {"-k", "-bpk", "-t"};
    const opt_set flags = {"-u", "-s", "--no-splice", "-f", "-v"};
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
                opt_map.contains("-f"), opt_map.contains("-v"),
                std::move(input), std::move(first_out), std::move(second_out));
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
    std::cerr << "  -v               output sizes of Bloom Filters" << std::endl;
    // clang-format on
    return 1;
}

std::string ComputeArgs::fasta_header() const {
    std::string mode = _unidirectional ? "unidirectional" : "bidirectional";
    std::string splice = _no_splice ? "false" : "true";
    std::stringstream ss;
    ss << "approximate masked superstring dataset='" << _dataset << "' k=" << _k
       << " bits-per-kmer=" << _bpk << " mode=" << mode << " splice=" << splice;
    return ss.str();
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
    std::string splice = _no_splice ? "false" : "true";
    std::stringstream ss;
    ss << "exact masked superstring dataset='" << _dataset << "' k=" << _k
       << " mode=" << mode << " splice=" << splice;
    return ss.str();
}

std::optional<CompareArgs> CompareArgs::from_cmdline(int argc,
                                                     std::string *argv) {
    const opt_set opts = {"-k"};
    const opt_set flags = {"-u"};
    std::unordered_map<std::string, std::string> opt_map = {};
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
        return CompareArgs(std::stoul(opt_map.at("-k")), opt_map.contains("-u"),
                           std::move(input), std::move(output));
    } catch (...) {
        return std::nullopt;
    }
}

int CompareArgs::usage() {
    // clang-format off
    std::cerr << "Usage: streaming-masked-superstrings compare [options] <approximate-output> <exact-output>" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -k <int>         kmer size" << std::endl;
    std::cerr << "  -u               treat kmer and its reverse complement as distinct" << std::endl;
    // clang-format on
    return 1;
}
