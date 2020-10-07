#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <bitset.h>
#include <exception>
#include <set>
#include <utility>
#include <map>

#define str second

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
using vb = std::vector<base::dynamic_bitset<>>;
using vusi = std::vector<std::unordered_set<int>>; // Adjacency list for indexed nodes

using umsb = um<std::string, base::dynamic_bitset<>>;

struct bimap_str_int {
    vs int_to_str;
    umsi str_to_int;

    int size() const { return int_to_str.size(); }

    bool has(const std::string &key) const { return str_to_int.find(key) != str_to_int.end(); };

    int index_of(const std::string &key) const { return str_to_int.at(key); };

    std::string key_at(const int index) const { return int_to_str[index]; };
};

using entities_bimap = bimap_str_int;

#endif // ! TYPES_HPP



