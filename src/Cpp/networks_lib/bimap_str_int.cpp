#include "bimap_str_int.hpp"

// Creates a bimap of string to int and viceversa.
// The int index assigned to each string corresponds to the lexicographic order.
// Creates a bimap of the elements in the selected column.
// Column index starts counting at 0
// The columns of the file must be separated by a tab ('\t')
bimap_str_int createBimap(std::string_view path_file, bool has_header, int column_index, int total_num_columns) {
    vs index_to_entities = createIntToStr(path_file, has_header, column_index, total_num_columns);
    umsi entities_to_index = createStrToInt(index_to_entities);
    return {index_to_entities, entities_to_index};
}

umsi createStrToInt(const vs &index_to_entities) {
    umsi entities_to_index;
    for (int I = 0; I < index_to_entities.size(); I++) {
        entities_to_index.emplace(index_to_entities[I], I);
    }
    return entities_to_index;
}

// Creates the bimap from a vector of elements
// Removes the duplicate elements in the vector
// Sorts the elements to assign the indexes
bimap_str_int createBimap(const vs &index_to_entities) {
    std::set<std::string> s(index_to_entities.begin(), index_to_entities.end());
    vs sortedAndUniqueVector(s.begin(), s.end());
    return {sortedAndUniqueVector, createStrToInt(sortedAndUniqueVector)};
}

// The input file has one identifier per row in the selected column.
// The index of the selected column starts counting at 0.
vs createIntToStr(std::string_view path_file, bool has_header, int selected_column, int total_num_columns) {
    std::ifstream map_file(path_file.data());
    std::string entity, leftover;
    uss temp_set;
    vs index_to_entities;
    int current_column;

    if (!map_file.is_open()) {
        std::string message = "Could not open file ";
        message += path_file;
        message += " at ";
        message += __FUNCTION__;
        throw std::runtime_error(message);
    }

    if (selected_column >= total_num_columns || selected_column < 0) {
        throw std::runtime_error("Invalid column index");
    }

    if (has_header) {
        getline(map_file, leftover);  // Skip header line
    }

    current_column = 0;
    while (map_file.peek() != EOF) {
        while (current_column != selected_column) {
            std::getline(map_file, leftover, '\t');
            current_column++;
        }
        if (selected_column == total_num_columns - 1) {
            std::getline(map_file, entity, '\n');
        } else {
            std::getline(map_file, entity, '\t');
            std::getline(map_file, leftover, '\n');
        }
        temp_set.insert(entity);
        current_column = 0;
    }
    index_to_entities = convert_uss_to_vs(temp_set);
    sort(index_to_entities.begin(), index_to_entities.end());

    return index_to_entities;
}