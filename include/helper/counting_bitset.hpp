#ifndef COUNTING_BITSET_HPP
#define COUNTING_BITSET_HPP

#include <algorithm>
#include <memory>

template <std::size_t BPC>
class CountingBitset {
  public:
    using inner_t = std::uint32_t;

    CountingBitset() : _size(0) {}
    CountingBitset(std::size_t size) : _size(size) {
        std::size_t inner_size = (size + cells_per_inner - 1) / cells_per_inner;
        data = std::make_unique<inner_t[]>(inner_size);
        std::fill(data.get(), data.get() + inner_size, 0);
    }
    void set(std::size_t ind, inner_t count) {
        if (count > max_count) {
            count = max_count;
        }
        reset(ind);
        std::size_t cell = get_index(ind);
        std::size_t offset = get_offset(ind);

        data[cell] |= (count << offset);
    }
    void reset(std::size_t ind) {
        std::size_t cell = get_index(ind);
        std::size_t offset = get_offset(ind);
        inner_t mask = ~(cell_mask << offset);

        data[cell] &= mask;
    }
    inner_t get(std::size_t ind) const {
        std::size_t cell = get_index(ind);
        std::size_t offset = get_offset(ind);
        return (data[cell] >> offset) & cell_mask;
    }
    bool test(std::size_t ind) const { return get(ind) > 0; }
    bool is_stuck(std::size_t ind) const { return get(ind) == max_count; }
    void increment(std::size_t ind) {
        if (get(ind) < max_count) {
            set(ind, get(ind) + 1);
        }
    }
    void decrement(std::size_t ind) {
        if (get(ind) > 0) {
            set(ind, get(ind) - 1);
        }
    }
    std::size_t size() const { return _size; }

  private:
    static constexpr std::size_t inner_size = sizeof(inner_t) * 8;
    static constexpr std::size_t max_count = (1 << BPC) - 1;
    static constexpr std::size_t cells_per_inner = inner_size / BPC;
    static constexpr inner_t cell_mask = (1 << BPC) - 1;

    static std::size_t get_index(std::size_t ind) {
        return ind / cells_per_inner;
    }
    static std::size_t get_offset(std::size_t ind) {
        return (ind % cells_per_inner) * BPC;
    }
    std::size_t _size;
    std::unique_ptr<inner_t[]> data;
};

#endif
