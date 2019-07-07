#ifndef BIMAP_HPP_
#define BIMAP_HPP_

#include <fstream>
#include <bitset>

#include "types.hpp"
#include "utility.hpp"

struct bimap {
	vs entities;
	umsi indexes;
};

bimap createBimap(std::string_view list_file_path, bool hasHeader = false);
bimap createBimap(const vs& index_to_entities);

vs createIndexToEntities(const std::string& list_file_path, bool hasHeader);
umsi createEntitiesToIndex(const vs& index_to_entities);

#endif // !BIMAP_HPP_
