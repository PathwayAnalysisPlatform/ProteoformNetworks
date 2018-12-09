#ifndef PATHWAY_ENTITY_H
#define PATHWAY_ENTITY_H

#include <algorithm>
#include <fstream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "proteoform.hpp"

namespace pathway {

enum entities { GENES,
                PROTEINS,
                PROTEOFORMS};

struct entities_bimap {
   std::vector<std::string> index_to_entities;
   std::unordered_map<std::string, int> entities_to_index;
};

std::vector<std::string> convert(const std::unordered_set<std::string>& a_set);

std::vector<std::string> createIndexToEntities(const std::string& entities_file_path);

std::unordered_map<std::string, int> createEntitiesToIndex(const std::vector<std::string>& index_to_entities);

entities_bimap readEntities(std::string_view entities_file_path);

}  // namespace pathway

#endif /* PATHWAY_ENTITY_H */