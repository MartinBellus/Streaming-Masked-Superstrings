#include "helper/bitset.hpp"

std::size_t align_up(std::size_t size, std::size_t align) {
    return (size + align - 1) & ~(align - 1);
}

DynamicBitset::DynamicBitset(std::size_t size) : _size(size) {
    std::size_t bit_size = align_up(size, inner_size);
    std::size_t inner_count = (bit_size + inner_size - 1) / inner_size;
    data = std::make_unique<inner_t[]>(inner_count);
    std::fill(data.get(), data.get() + inner_count, 0);
}

void DynamicBitset::set(std::size_t ind) {
    if (ind >= _size) {
        return;
    }
    data[ind / inner_size] |= (1 << (ind % inner_size));
}

void DynamicBitset::reset(std::size_t ind) {
    if (ind >= _size) {
        return;
    }
    data[ind / inner_size] &= ~(1 << (ind % inner_size));
}

bool DynamicBitset::test(std::size_t ind) const {
    if (ind >= _size) {
        return false;
    }
    return data[ind / inner_size] & (1 << (ind % inner_size));
}
