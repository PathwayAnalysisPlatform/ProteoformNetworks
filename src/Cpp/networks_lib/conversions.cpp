#include "conversions.hpp"

// The result vector is sorted lexicographycally
vs convert_uss_to_vs(const uss &a_set) {
    vs result;
    result.assign(a_set.begin(), a_set.end());
    return result;
}

// Reads the mapping from one string to another type of string objects.
// The mapping file must follow the convention that there are two columns separated by a tab
// The second column can contain multiple values separated by a blank space
// The first line of the file can have the header with the name of the columns.
entity_mapping readMapping(std::string_view path_file_mapping, bool has_header_row) {
    entity_mapping mapping;
    std::ifstream file_mapping(path_file_mapping.data());
    std::string source, destination, destinations, leftover;

    if (!file_mapping.is_open()) {
        std::string message = "Error reading file at: ";
        message += __FUNCTION__;
        throw std::runtime_error(message);
    }

    if (has_header_row)
        getline(file_mapping, leftover);    // Discard the header line.
    while (getline(file_mapping, source, '\t')) {
        while (getline(file_mapping, destinations)) {   // get the whole second column until the end of the line
            std::stringstream ss(destinations);
            while (getline(ss, destination, ' ')) {
                mapping.first_to_second.emplace(source, destination);
                mapping.second_to_first.emplace(destination, source);
            }
        }
    }

    std::cerr << "Mapping loaded: " << mapping.first_to_second.size() << " = "
              << mapping.second_to_first.size() << "\n";
    return mapping;
}

