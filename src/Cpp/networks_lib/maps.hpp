#ifndef PROTEOFORMNETWORKS_MAPS_HPP
#define PROTEOFORMNETWORKS_MAPS_HPP

#include <set>
#include "types.hpp"

template<class K, class V>
bool hasKey(const std::map<K, V> &m, K key) {
    return m.find(key) != m.end();
}

template<class K, class V>
bool hasKey(const um<K, V> &m, K key) {
    return m.find(key) != m.end();
}

template<class K, class V>
bool keyHasValue(const um<K, V> &m, K key, V value) {
    if (!hasKey(m, key))
        return false;
    return m[key] == value;
}

template<class K, class V>
bool hasKey(const umm<K, V> &m, K key) {
    return m.find(key) != m.end();
}

template<class K, class V>
bool keyHasValue(const umm<K, V> &m, K key, V value) {
    auto range = m.equal_range(key);
    for (auto it = range.first; it != range.second; it++) {
        if (it->second == value)
            return true;
    }
    return false;
}

template<class K, class V>
bool hasValue(const umm<K, V> &m, V value) {
    for (auto it = m.begin(); it != m.end(); it++)
        if (it->second == value)
            return true;
    return false;
}

template<class K, class V>
bool hasValue(const um<K, V> &m, V value) {
    for (auto it = m.begin(); it != m.end(); it++)
        if (it->second == value)
            return true;
    return false;
}

template<class V>
bool hasValue(const us<V> &s, V value) {
    return s.find(value) != s.end();
}

template<class V>
bool hasValue(const std::set<V> &s, V value) {
    return s.find(value) != s.end();
}

template<class V>
bool hasValue(const std::vector<V> &v, V value) {
    for (const auto &element : v) {
        if (element == value)
            return true;
    }
    return false;
}

template<class T, class V>
std::unordered_set<T> getKeys(const std::unordered_multimap<T, V> &m) {
    std::unordered_set<T> keys;
    for (auto it = m.begin(); it != m.end(); it++)
        keys.insert(it->first);
    return keys;
}

template<class T, class V>
std::unordered_set<T> getKeys(const std::map<T, V> &m) {
    std::unordered_set<T> keys;
    for (auto it = m.begin(); it != m.end(); it++)
        keys.insert(it->first);
    return keys;
}

template<class T, class V>
std::unordered_set<T> getValues(const std::unordered_multimap<T, V> &m) {
    std::unordered_set<V> values;
    for (auto it = m.begin(); it != m.end(); it++)
        values.insert(it->second);
    return values;
}

#endif //PROTEOFORMNETWORKS_MAPS_HPP
