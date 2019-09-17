#include "networks.hpp"

// Read entity interactions in Reactome from a PathwayMatcher edges file
// For each interacting entities A and B, adds both the edges A --> B and B --> A
ummii loadInteractionNetwork(std::string_view path_file_interactions,
                             const bimap_str_int &entities,
                             bool has_header_row) {

    ummii interactions;
    std::ifstream file_interactions(path_file_interactions.data());
    std::string e1, e2, other_fields;

    if (!file_interactions.is_open()) {
        std::string message = "Error reading file at: ";
        message += __FUNCTION__;
        throw std::runtime_error(message);
    }

    if (has_header_row)
        getline(file_interactions, other_fields);    // Discard the header line.
    while (getline(file_interactions, e1, '\t')) {
        getline(file_interactions, e2, '\t');   // get second entity
        getline(file_interactions, other_fields);// read other columns
        if (entities.str_to_int.find(e1) != entities.str_to_int.end()
            && entities.str_to_int.find(e2) != entities.str_to_int.end()) {
            int index_e1 = entities.str_to_int.at(e1);
            int index_e2 = entities.str_to_int.at(e2);
            interactions.emplace(index_e1, index_e2);
            interactions.emplace(index_e2, index_e1);
        } else {
            if(entities.str_to_int.find(e1) == entities.str_to_int.end()){
                std::cerr << "Not found entity: " << e1 << "\n";
            }
            if(entities.str_to_int.find(e2) == entities.str_to_int.end()){
                std::cerr << "Not found entity: " << e2 << "\n";
            }
        }
    }

    return interactions;
}