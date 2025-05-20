#ifndef ARGS_HPP
#define ARGS_HPP

#include <optional>
#include <string>

class ComputeArgs {
  public:
    static std::optional<ComputeArgs> from_cmdline(int argc, std::string *argv);
    static int usage();
    std::size_t k() const { return _k; }
    std::size_t bits_per_element() const { return _bpk; }
    const std::string &dataset() const { return _dataset; }
    const std::string &output() const { return _output; }

    std::string fasta_header() const;

  private:
    ComputeArgs(std::size_t k, std::size_t bpk, std::string &&dataset,
                std::string &&output)
        : _k(k), _bpk(bpk), _dataset(std::move(dataset)),
          _output(std::move(output)) {}
    std::size_t _k;
    std::size_t _bpk;
    std::string _dataset;
    std::string _output;
};

class ExactArgs {
  public:
    static std::optional<ExactArgs> from_cmdline(int argc, std::string *argv);
    static int usage();

    std::size_t k() const { return _k; }
    const std::string &dataset() const { return _dataset; }
    const std::string &output() const { return _output; }

    std::string fasta_header() const;

  private:
    ExactArgs(std::size_t k, std::string &&dataset, std::string &&output)
        : _k(k), _dataset(std::move(dataset)), _output(std::move(output)) {}
    std::size_t _k;
    std::string _dataset;
    std::string _output;
};

class CompareArgs {
  public:
    static std::optional<CompareArgs> from_cmdline(int argc, std::string *argv);
    static int usage();

    const std::string &output() const { return _output; }
    const std::string &golden() const { return _golden; }

  private:
    CompareArgs(std::string &&output, std::string &&golden)
        : _output(std::move(output)), _golden(std::move(golden)) {}
    std::string _output;
    std::string _golden;
};

#endif
