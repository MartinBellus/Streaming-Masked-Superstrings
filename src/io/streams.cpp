#include "io/streams.hpp"

using namespace io;

char input_stream::peek() { return stream.peek(); }

bool input_stream::eof() const { return stream.eof(); }

bool input_stream::getline(std::string &s) {
    return std::getline(stream, s).good();
}

bool input_stream::get(char &c) { return stream.get(c).good(); }

void input_stream::ignore() { stream.ignore(); }

void output_stream::write(char c) { stream.put(c); }
