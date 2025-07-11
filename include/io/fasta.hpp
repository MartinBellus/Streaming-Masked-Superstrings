#ifndef FASTA_HPP
#define FASTA_HPP

#include "helper/kmer.hpp"
#include "io/streams.hpp"
#include <string>

namespace io {

class FastaReader {
  private:
    static constexpr char CommentChar = '>';
    static constexpr std::size_t MaxHeaderSize = 80;

  public:
    FastaReader(input_stream &&stream) : stream(std::move(stream)) {}
    FastaReader(const std::string &path) : stream(path) {}
    bool next_sequence() {
        if (stream.eof() || stream.peek() != CommentChar)
            return false;
        stream.ignore();
        return stream.getline(header);
    }
    bool next_nucleotide(char &next) {
        skip_ws();
        if (stream.eof() || stream.peek() == CommentChar) {
            return false;
        }
        return stream.get(next);
    }
    const std::string &get_header() const { return header; }
    void reset() {
        stream.reset();
        header.clear();
    }

  private:
    void skip_ws() {
        char c;
        while ((c = stream.peek()) && std::isspace(c)) {
            stream.ignore();
        }
    }
    input_stream stream;
    std::string header;
};

enum {
    NOT_PRESENT = 0,
    PRESENT = 1,
};

class KmerWriter {
  public:
    KmerWriter(output_stream &&stream, std::size_t K, bool splice)
        : stream(std::move(stream)), kmer(K), last_one(K), splice(splice) {}
    KmerWriter(const std::string &path, std::size_t K, bool splice)
        : stream(path), kmer(K), last_one(K), splice(splice) {}
    void write_header(const std::string &header);
    void add_nucleotide(char c) { kmer.roll(c); }
    void print_nucleotide(int present);
    void flush();

  private:
    output_stream stream;
    Kmer kmer;
    std::size_t last_one;
    bool splice;
};

} // namespace io

#endif
