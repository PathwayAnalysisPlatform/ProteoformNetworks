#include <cstring>
#include "bimap_str_int.hpp"

using namespace std;

umsi createStrToInt(const vs& index_to_entities) {
	umsi entities_to_index;
	for (int I = 0; I < index_to_entities.size(); I++) {
		entities_to_index.emplace(index_to_entities[I], I);
	}
	return entities_to_index;
}

// Creates a bimap of string to int and viceversa.
// The int index assigned to each string corresponds to the lexicographic order.
bimap_str_int createBimap(string_view list_file_path, bool has_header) {
	vs index_to_entities = createIntToStr(list_file_path, has_header);
	umsi entities_to_index = createStrToInt(index_to_entities);
	return { index_to_entities, entities_to_index };
}

bimap_str_int createBimap(const vs& index_to_entities) {
	umsi entities_to_index = createStrToInt(index_to_entities);
	return { index_to_entities, entities_to_index };
}

// The input file is a list of identifiers. One in each row
vs createIntToStr(string_view list_file_path, bool has_header) {
	ifstream map_file(list_file_path.data());
	string entity, leftover;
	uss temp_set;
	vs index_to_entities;

	if (!map_file.is_open()) {
		std::string message = "Could not open file list_file_path ";
        std::string function = __FUNCTION__;
        throw runtime_error(message + function);
	}

	if (has_header) {
		getline(map_file, leftover);  // Skip header line
	}
	while (map_file.peek() != EOF) {
		getline(map_file, entity);
		temp_set.insert(entity);
	}
	index_to_entities = convert_uss_to_vs(temp_set);
    sort(index_to_entities.begin(), index_to_entities.end());

	return index_to_entities;
}