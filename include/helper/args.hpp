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
    bool unidirectional() const { return _unidirectional; }
    bool splice() const { return !_no_splice; }
    bool second_phase() const { return !_skip_second_phase; }
    bool verbose() const { return _verbose; }
    const std::string &dataset() const { return _dataset; }
    const std::string &first_phase_output() const { return _first_out; }
    const std::string &second_phase_output() const { return _second_out; }

    std::string fasta_header() const;

  private:
    ComputeArgs(std::size_t k, std::size_t bpk, bool unidirectional,
                bool splice, bool skip_second, bool verbose,
                std::string &&dataset, std::string &&first_out,
                std::string &&second_out)
        : _k(k), _bpk(bpk), _unidirectional(unidirectional), _no_splice(splice),
          _skip_second_phase(skip_second), _verbose(verbose),
          _dataset(std::move(dataset)), _first_out(std::move(first_out)),
          _second_out(std::move(second_out)) {}
    std::size_t _k;
    std::size_t _bpk;
    bool _unidirectional;
    bool _no_splice;
    bool _skip_second_phase;
    bool _verbose;
    std::string _dataset;
    std::string _first_out;
    std::string _second_out;
};

class ExactArgs {
  public:
    static std::optional<ExactArgs> from_cmdline(int argc, std::string *argv);
    static int usage();

    std::size_t k() const { return _k; }
    bool unidirectional() const { return _unidirectional; }
    bool splice() const { return !_no_splice; }
    const std::string &dataset() const { return _dataset; }
    const std::string &output() const { return _output; }

    std::string fasta_header() const;

  private:
    ExactArgs(std::size_t k, bool unidirectional, bool splice,
              std::string &&dataset, std::string &&output)
        : _k(k), _unidirectional(unidirectional), _no_splice(splice),
          _dataset(std::move(dataset)), _output(std::move(output)) {}
    std::size_t _k;
    bool _unidirectional;
    bool _no_splice;
    std::string _dataset;
    std::string _output;
};

class CompareArgs {
  public:
    static std::optional<CompareArgs> from_cmdline(int argc, std::string *argv);
    static int usage();

    std::size_t k() const { return _k; }
    bool unidirectional() const { return _unidirectional; }
    const std::string &output() const { return _output; }
    const std::string &golden() const { return _golden; }

  private:
    CompareArgs(std::size_t k, bool unidirectional, std::string &&output,
                std::string &&golden)
        : _k(k), _unidirectional(unidirectional), _output(std::move(output)),
          _golden(std::move(golden)) {}
    std::size_t _k;
    bool _unidirectional;
    std::string _output;
    std::string _golden;
};

#endif
