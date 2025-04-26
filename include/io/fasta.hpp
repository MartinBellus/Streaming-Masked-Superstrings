#ifndef FASTA_HPP
#define FASTA_HPP

#include <string>

template <typename T>
class FastaReader {
  private:
    static constexpr char CommentChar = '>';
    static constexpr std::size_t MaxHeaderSize = 80;

  public:
    FastaReader(T &&stream) : stream(std::forward<T>(stream)) {}
    bool next_sequence() {
        if (stream.eof() || stream.peek() != CommentChar)
            return false;
        stream.ignore();
        return std::getline(stream, header).good();
    }
    bool next_nucleotide(char &next) {
        skip_ws();
        if (stream.eof() || stream.peek() == CommentChar) {
            return false;
        }
        return stream.get(next).good();
    }
    const std::string &get_header() const { return header; }

  private:
    void skip_ws() {
        char c;
        while ((c = stream.peek()) && std::isspace(c)) {
            stream.ignore();
        }
    }
    T stream;
    std::string header;
};

#endif
