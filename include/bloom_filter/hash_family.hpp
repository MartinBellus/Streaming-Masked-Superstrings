#ifndef HASH_FAMILY_HPP
#define HASH_FAMILY_HPP

#include <memory>

class hash_family {
  public:
    using hash_t = std::uint64_t;
    std::span<const hash_t> hash(this auto &&self, const std::string &s) {
        return self.hash_impl(s);
    }
    hash_family(std::size_t nhashes) : nhashes(nhashes) {
        buffer = std::make_unique<hash_t[]>(nhashes);
    }

  protected:
    std::unique_ptr<hash_t[]> buffer;
    std::size_t nhashes;
};

class rolling_hash_family : public hash_family {
  public:
    using hash_t = std::uint64_t;

    rolling_hash_family(std::size_t nhashes) : hash_family(nhashes) {}
    void roll(this auto &&self, char c) { self.roll_impl(c); }
    void init(this auto &&self, const std::string &s) { self.init_impl(s); }
    std::span<const hash_t> hash_impl(this auto &&self, const std::string &s) {
        self.init(s);
        return self.get_hashes();
    }
    std::span<const hash_t> get_hashes() const {
        return std::span(buffer.get(), nhashes);
    }
};

template <class T>
concept HashFamily = std::derived_from<T, hash_family> && requires(T t) {
    {
        t.hash_impl(std::declval<std::string &>())
    } -> std::same_as<std::span<const typename T::hash_t>>;
};

template <class T>
concept RollingHashFamily =
        std::derived_from<T, rolling_hash_family> && requires(T t) {
            { t.roll_impl(std::declval<char>()) } -> std::same_as<void>;
            {
                t.init_impl(std::declval<std::string &>())
            } -> std::same_as<void>;
            {
                t.get_hashes()
            } -> std::same_as<std::span<const typename T::hash_t>>;
        };

#endif
