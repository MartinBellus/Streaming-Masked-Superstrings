#ifndef BITSET_HPP
#define BITSET_HPP

#include <cstdint>
#include <memory>

class DynamicBitset {
  private:
    using inner_t = std::uint32_t;
    static constexpr std::size_t inner_size = sizeof(inner_t) * 8;

  public:
    DynamicBitset() : _size(0) {}
    DynamicBitset(std::size_t size);
    void set(std::size_t ind);
    void reset(std::size_t ind);
    bool test(std::size_t ind) const;
    std::size_t size() const { return _size; }

  private:
    std::size_t _size;
    std::unique_ptr<inner_t[]> data;
};

#endif
