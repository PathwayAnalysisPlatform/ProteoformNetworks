#ifndef PROTEOFORMNETWORKS_MAPS_HPP
#define PROTEOFORMNETWORKS_MAPS_HPP

#include <set>
#include "types.hpp"

bool hasKey(const umsi &m, std::string key);
bool hasKey(const ummii &m, int key);
bool keyHasValue(const ummii &m, int key, int value);

template<class T, class V>
bool hasKey(const um<T, V> &m, T key) {
    return m.find(key) != m.end();
}

template<class T>
bool hasKey(const us<T> &s, T key) {
    return s.find(key) != s.end();
}

#endif //PROTEOFORMNETWORKS_MAPS_HPP
