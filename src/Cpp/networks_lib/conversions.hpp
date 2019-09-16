#ifndef PROTEOFORMNETWORKS_CONVERSIONS_HPP
#define PROTEOFORMNETWORKS_CONVERSIONS_HPP

#include "types.hpp"

#include <algorithm>

// Convert unordered string set, to string vector
vs convert_uss_to_vs(const uss &a_set);

entity_mapping readMapping(std::string_view path_file_mapping);

#endif //PROTEOFORMNETWORKS_CONVERSIONS_HPP