#ifndef BIMAP_HPP_
#define BIMAP_HPP_

#include <fstream>

#include "types.hpp"
#include "utility.hpp"

struct bimap {
	vs entities;
	umsi indexes;
};

vs createIndexToEntities(const std::string& path_file, bool hasHeader = false);

umsi createEntitiesToIndex(const vs& index_to_entities);

bimap createBimap(std::string_view entities_file_path, bool hasHeader = false);
bimap createBimap(vs index_to_entities);

#endif // !BIMAP_HPP_
