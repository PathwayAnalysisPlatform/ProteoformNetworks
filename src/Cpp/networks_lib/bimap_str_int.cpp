#include "bimap_str_int.hpp"

std::string rtrim(std::string &s) {
    if (s.length() > 0) {
        auto it = s.end() - 1;
        while (isspace(*it)) {
            s.erase(it);
            it = s.end() - 1;
        }
    }
    return s;
}

vs convert_uss_to_vs(const uss &a_set) {
    vs result;
    result.assign(a_set.begin(), a_set.end());
    return result;
}

// Create from list, without header, list comes already sorted, there is only one column in the file.
Bimap_str_int::Bimap_str_int(std::string_view file_elements) {
    std::cout << "Reading vertices...\n";
    std::ifstream f;
    f.open(file_elements.data());

    if (!f.is_open()) {
        std::string message = "Cannot open vertices file at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    std::string element;
    int i = 0;
    while (f >> element) {
//        std::cout << "Adding element: " << element << std::endl;
        itos.push_back(element);
        stoi[element] = i;
        i++;
    }
    std::cout << "Complete.\n\n";
}

// Creates a bimap of string to int and viceversa.
// The int index assigned to each string corresponds to the lexicographic order.
// Creates a bimap of the elements in the selected column.
// Column index starts counting at 0
// The columns of the file must be separated by a tab ('\t')
Bimap_str_int::Bimap_str_int(std::string_view file_elements, bool has_header, int column_index, int total_num_columns) :
        itos{createIntToStr(file_elements, has_header, column_index, total_num_columns)},
        stoi{createStrToInt(itos)} {
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
Bimap_str_int::Bimap_str_int(const vs &index_to_entities) {
    std::set<std::string> s(index_to_entities.begin(), index_to_entities.end());
    vs sortedAndUniqueVector(s.begin(), s.end());
    itos = sortedAndUniqueVector;
    stoi = createStrToInt(sortedAndUniqueVector);
}

int Bimap_str_int::index(const std::string &key) const {
    if(stoi.find(key) == stoi.end())
    {
        std::cout << "Key not found: " << key << std::endl;
        return -1;
    }
    return stoi.at(key);
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
        temp_set.insert(rtrim(entity));
        current_column = 0;
    }
    index_to_entities = convert_uss_to_vs(temp_set);
    sort(index_to_entities.begin(), index_to_entities.end());

    return index_to_entities;
}


