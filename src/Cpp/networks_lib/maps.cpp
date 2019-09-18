#include "maps.hpp"

bool hasKey(const umsi &m, std::string key) {
    return m.find(key) != m.end();
}

bool hasKey(const ummii &m, int key) {
    return m.find(key) != m.end();
}

bool keyHasValue(const ummii &m, int key, int value) {
    auto range = m.equal_range(key);
    for(auto it = range.first; it != range.second; it++){
        if(it->second == value)
            return true;
    }
    return false;
}
