#ifndef PROTEOFORMNETWORKS_CONVERSIONS_HPP
#define PROTEOFORMNETWORKS_CONVERSIONS_HPP

#include "types.hpp"

#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <sstream>
#include <iostream>

// Convert unordered string set, to string vector
vs convert_uss_to_vs(const uss &a_set);

entity_mapping readMapping(std::string_view path_file_mapping, bool has_header_row = true);

#endif //PROTEOFORMNETWORKS_CONVERSIONS_HPP