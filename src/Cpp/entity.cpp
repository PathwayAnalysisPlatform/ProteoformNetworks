#include "entity.hpp"

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

std::vector<std::string> createIndexToEntities(const std::string& path_file_mapping) {
   std::ifstream map_file(path_file_mapping.data());
   std::string entity, leftover;
   std::unordered_set<std::string> temp_set;
   std::vector<std::string> index_to_entities;

   if (!map_file.is_open()) {
      throw std::runtime_error("Could not open file " + path_file_mapping);
   }

   getline(map_file, leftover);  // Skip header line
   while (map_file.peek() != EOF) {
      if (map_file.peek() == '[' || map_file.peek() == '\"') {
         if (map_file.peek() == '\"')  // Read initial "
            map_file.get();
         map_file.get();  // Read initial " or [
         getline(map_file, entity, ']');
      } else {
         getline(map_file, entity, ',');  // Read entity
      }
      getline(map_file, leftover);  // Read rest of line
      temp_set.insert(entity);
   }
   index_to_entities = convert(temp_set);

   return index_to_entities;
}

entities_bimap readEntities(std::string_view path_file_mapping) {
   auto index_to_entities = createIndexToEntities(path_file_mapping.data());
   auto entities_to_index = createEntitiesToIndex(index_to_entities);
   return {index_to_entities, entities_to_index};
}

}  // namespace pathway
