#ifndef BIMAP_HPP_
#define BIMAP_HPP_

#include <fstream>
#include <bitset>
#include <string_view>
#include <cstring>
#include <set>
#include "conversions.hpp"

// Create bimap from element to index in a unique sorted array.
// Creates a bimap of the elements in the selected column.
// Column index starts counting at 0
bimap_str_int
createBimap(std::string_view path_file, bool has_header = true, int column_index = 0, int total_num_columns = 1);

bimap_str_int createBimap(const vs &index_to_entities);

vs createIntToStr(std::string_view path_file, bool has_header = true, int column_index = 0, int total_num_columns = 1);

umsi createStrToInt(const vs &index_to_entities);

#endif // !BIMAP_HPP_
