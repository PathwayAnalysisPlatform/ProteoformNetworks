#ifndef PROTEOFORMNETWORKS_CONVERSIONS_HPP
#define PROTEOFORMNETWORKS_CONVERSIONS_HPP

#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <sstream>
#include <iostream>

#include "types.hpp"
#include "overlap_types.hpp"

// Convert unordered string set, to string vector
vs convert_uss_to_vs(const uss &a_set);

entity_mapping
readMapping(std::string_view path_file_mapping, bool has_header_row = true, bool has_additional_columns = false);

#endif //PROTEOFORMNETWORKS_CONVERSIONS_HPP