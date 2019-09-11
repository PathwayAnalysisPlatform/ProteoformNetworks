#ifndef BIMAP_HPP_
#define BIMAP_HPP_

#include <fstream>
#include <bitset>
#include <string_view>

#include "../others/types.hpp"
#include "conversions.hpp"

struct bimap_str_int {
	vs int_to_str;
	umsi str_to_int;
};

bimap_str_int createBimap(std::string_view list_file_path, bool has_header = true);
bimap_str_int createBimap(const vs& index_to_entities);

vs createIntToStr(std::string_view list_file_path, bool has_header = true);
umsi createStrToInt(const vs& index_to_entities);

#endif // !BIMAP_HPP_
