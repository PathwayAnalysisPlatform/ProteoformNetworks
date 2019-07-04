#include "bimap.hpp"

umsi createEntitiesToIndex(const vs& index_to_entities) {
	umsi entities_to_index;
	for (int I = 0; I < index_to_entities.size(); I++) {
		entities_to_index.emplace(index_to_entities[I], I);
	}
	return entities_to_index;
}

// The input file is a list of identifiers. One in each row
vs createIndexToEntities(const std::string& path_file, bool hasHeader) {
	std::ifstream map_file(path_file.data());
	std::string entity, leftover;
	uss temp_set;
	vs index_to_entities;

	if (!map_file.is_open()) {
		throw std::runtime_error("Could not open file " + path_file);
	}

	if (hasHeader) {
		getline(map_file, leftover);  // Skip header line
	}
	while (map_file.peek() != EOF) {
		getline(map_file, entity);
		temp_set.insert(entity);
	}
	index_to_entities = convert_uss_to_vs(temp_set);

	return index_to_entities;
}

bimap createBimap(std::string_view path_file_mapping, bool hasHeader) {
	vs index_to_entities = createIndexToEntities(path_file_mapping.data(), hasHeader);
	umsi entities_to_index = createEntitiesToIndex(index_to_entities);
	return { index_to_entities, entities_to_index };
}

bimap createBimap(vs index_to_entities) {
	umsi entities_to_index = createEntitiesToIndex(index_to_entities);
	return { index_to_entities, entities_to_index };
}