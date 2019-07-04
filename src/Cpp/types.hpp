#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using msi = std::map<std::string, int>;

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

using vs = std::vector<std::string>;

template<typename K>
using bm = std::pair<std::vector<K>, um<K, int>>;
using bms = bm<std::string>;	// TODO: Replace bimap for this

#endif // ! TYPES_HPP



