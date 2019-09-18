#ifndef PROTEOFORMNETWORKS_MAPS_HPP
#define PROTEOFORMNETWORKS_MAPS_HPP

#include <set>
#include "types.hpp"

bool hasKey(const umsi &m, std::string key);
bool hasKey(const ummii &m, int key);

bool keyHasValue(const ummii &m, int key, int value);

#endif //PROTEOFORMNETWORKS_MAPS_HPP
