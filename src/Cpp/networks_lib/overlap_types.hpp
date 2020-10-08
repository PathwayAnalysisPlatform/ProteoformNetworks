#ifndef PROTEOFORMNETWORKS_OVERLAP_TYPES_HPP
#define PROTEOFORMNETWORKS_OVERLAP_TYPES_HPP

#include <types.hpp>

// A hash function used to hash a pair of any kind
struct hash_pair {
    template<class T1, class T2>
    size_t operator()(const std::pair<T1, T2> &p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

template<typename T>
using pair_map = std::unordered_map<std::pair<int, int>, T, hash_pair>;

#endif //PROTEOFORMNETWORKS_OVERLAP_TYPES_HPP
