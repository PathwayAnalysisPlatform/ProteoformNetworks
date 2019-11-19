#include "conversions.hpp"

// The result vector is sorted lexicographycally
vs convert_uss_to_vs(const uss &a_set) {
    vs result;
    result.assign(a_set.begin(), a_set.end());
    return result;
}

// Reads the mapping from one string to another type of string objects.
// THere are two possible formats, both with tab separator for the columns:
//  1) Two columns, first one has one id, the second column contains one or more values separated by a blank space
//  2) Multiple columns, the first two one id, and the second only one other id.
//  The first line of the file can have the header with the name of the columns.
// There should be no extra spaces after the last string of the row
bidirectional_mapping
readMapping(std::string_view path_file_mapping, bool has_header_row, bool has_additional_columns) {
    bidirectional_mapping mapping;
    std::ifstream file_mapping(path_file_mapping.data());
    std::string source, destination, destinations, leftover;

    if (!file_mapping.is_open()) {
        std::string message = "Error reading file ";
        message += path_file_mapping;
        message += " at: ";
        message += __FUNCTION__;
        throw std::runtime_error(message);
    }

    if (has_header_row)
        std::getline(file_mapping, leftover);    // Discard the header line.

    while (std::getline(file_mapping, source, '\t')) {

        if (has_additional_columns) {
            std::getline(file_mapping, destination, '\t');
            std::getline(file_mapping, leftover, '\n');
            mapping.first_to_second.emplace(source, destination);
            mapping.second_to_first.emplace(destination, source);
        } else {
            std::getline(file_mapping, destinations, '\n');   // get the whole second column until the end of the line
            std::stringstream ss(destinations);
            while (std::getline(ss, destination, ' ')) {
                mapping.first_to_second.emplace(source, destination);
                mapping.second_to_first.emplace(destination, source);
            }
        }
    }

    std::cerr << "Mapping loaded: " << mapping.first_to_second.size() << " = "
              << mapping.second_to_first.size() << "\n";
    return mapping;
}

