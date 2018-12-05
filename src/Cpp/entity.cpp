#include "entity.h"

namespace pathway {

std::vector<std::string> convert(const std::unordered_set<std::string>& a_set) {
   std::vector<std::string> result;
   result.assign(a_set.begin(), a_set.end());
   sort(result.begin(), result.end());
   return result;
}

std::unordered_map<std::string, int> createEntitiesToIndex(const std::vector<std::string>& index_to_entities) {
   std::unordered_map<std::string, int> entities_to_index;
   for (int I = 0; I < index_to_entities.size(); I++) {
      entities_to_index.emplace(index_to_entities[I], I);
   }
   return entities_to_index;
}

std::vector<std::string> createIndexToEntities(const std::string& entities_file_path) {
   std::ifstream entities_file(entities_file_path.data());
   std::string entity, line_leftover;
   std::unordered_set<std::string> temp_set;
   std::vector<std::string> index_to_entities;

   if (!entities_file.is_open()) {
      throw std::runtime_error("Could not open file " + entities_file_path);
   }

   while (getline(entities_file, entity, '\t')) {
      temp_set.insert(entity);
      getline(entities_file, line_leftover);
   }
   index_to_entities = convert(temp_set);

   return index_to_entities;
}

entities_bimap readEntities(std::string_view entities_file_path) {
   auto index_to_entities = createIndexToEntities(entities_file_path.data());
   auto entities_to_index = createEntitiesToIndex(index_to_entities);
   return {index_to_entities, entities_to_index};
}

}  // namespace pathway
