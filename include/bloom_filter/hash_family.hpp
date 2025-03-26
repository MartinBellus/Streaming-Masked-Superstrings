#ifndef HASH_FAMILY_HPP
#define HASH_FAMILY_HPP

#include <type_traits>

class hash_family_tag {};

template <class T>
concept HashFamily = std::is_base_of_v<hash_family_tag, T>;

#endif
