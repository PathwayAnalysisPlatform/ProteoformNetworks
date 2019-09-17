#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <bitset.h>

using msi = std::map<std::string, int>;
using msb = std::map<std::string, base::dynamic_bitset<>>;

template<typename K>
using us = std::unordered_set<K>;
using uss = us<std::string>;

template<typename K, typename V>
using um = std::unordered_map<K, V>;
using umss = um<std::string, std::string>;
using umsi = um<std::string, int>;

template<typename K, typename V>
using umm = std::unordered_multimap<K, V>;
using ummss = umm<std::string, std::string>;
using ummsi = std::unordered_multimap<std::string, int>;
using ummii = std::unordered_multimap<int, int>;

using vs = std::vector<std::string>;

using umsb = um<std::string, base::dynamic_bitset<>>;

struct entity_mapping {
    ummss first_to_second;
    ummss second_to_first;
};

struct trait_modules {
    msb traits_to_entities;
    msb entities_to_traits;
};

#endif // ! TYPES_HPP



