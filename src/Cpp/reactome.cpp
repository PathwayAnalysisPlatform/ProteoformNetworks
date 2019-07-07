#include "reactome.hpp"

// Read list of entities from a PathwayMatcher search file
vs loadReactomeEntitiesIndexToEntitites(std::string_view path_search_file) {
	vs result;
	std::ifstream file_search(path_search_file.data());
	std::string field, entity;
	uss temp_set;

	getline(file_search, field);					// Skip header line
	while (getline(file_search, entity, '\t')) {	// Read Entity (Gene | Protein | Proteoform)
		getline(file_search, field);				// Read rest of line
		temp_set.insert(entity);
	}
	return convert_uss_to_vs(temp_set);
}

// Create bimap from entities in a PathwayMatcher search file
bimap loadReactomeEntities(std::string_view path_search_file) {
	auto index_to_entities = loadReactomeEntitiesIndexToEntitites(path_search_file.data());
	return createBimap(index_to_entities);
}