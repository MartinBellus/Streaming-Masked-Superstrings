#ifndef STREAMS_HPP
#define STREAMS_HPP

#include <fstream>
#include <string>

namespace io {

class input_stream {
  public:
    input_stream(std::ifstream &&stream) : stream(std::move(stream)) {}
    input_stream(const std::string &path) : stream(path) {}
    bool is_open() const { return stream.is_open(); }
    char peek();
    bool eof() const;
    bool getline(std::string &);
    bool get(char &);
    void ignore();
    void reset();

  private:
    std::ifstream stream;
};

class output_stream {
  public:
    output_stream(std::ofstream &&stream) : stream(std::move(stream)) {}
    output_stream(const std::string &path) : stream(path) {}
    void write(char c);
    void write(const std::string &);
    bool is_open() const { return stream.is_open(); }

  private:
    std::ofstream stream;
};

} // namespace io
#endif
