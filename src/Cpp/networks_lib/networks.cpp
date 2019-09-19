#include "networks.hpp"
#include "overlap_types.hpp"

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

modules removeDisconnectedMembers(modules modules, ummii interactions) {
    // For each module:
        // Get module members

        // For each member, check if it has a neighbor in the members

        // Check if each member interacts with any other member of the module


    return modules;
}

// Reads groups from files. It creates the bimap structures for the groups and the members
// The file should follow the format convention: two colums, no extra spaces at the end of each row.
// Columns separated by a '\t' (tab) character.
load_modules_result loadModules(std::string_view path_file_modules, bool has_header) {
    modules result_modules;
    bimap_str_int groups;
    bimap_str_int members;

    std::ifstream file_modules(path_file_modules.data());
    if (!file_modules.is_open()) {
        std::string message = "Cannot open path_file_modules at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    // # Read the file once to get the groups and members lists
    std::string group, member, line;
    std::set<std::string> unique_groups, unique_members;
    if (has_header)
        std::getline(file_modules, line);
    while (std::getline(file_modules, group, '\t')) {
        getline(file_modules, member);
        unique_groups.insert(group);
        unique_members.insert(member);
    }
    groups = createBimap(vs(unique_groups.begin(), unique_groups.end()));
    members = createBimap(vs(unique_members.begin(), unique_members.end()));
    file_modules.close();

    // # Read the file to create the modules
    file_modules.open(path_file_modules.data());
    if (!file_modules.is_open()) {
        std::string message = "Cannot open path_file_modules at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    // ## Read modules
    // Initialize the modules to the correct sizes
    for (const auto &group : groups.int_to_str)
        result_modules.group_to_members.emplace(group, base::dynamic_bitset<>(groups.int_to_str.size()));
    for (const auto &member : members.int_to_str)
        result_modules.member_to_groups.emplace(member, base::dynamic_bitset<>(members.int_to_str.size()));

    // Set the members of each group and the owners of each member
    if (has_header)
        std::getline(file_modules, line);
    while (std::getline(file_modules, group, '\t')) {
        std::getline(file_modules, member);
        result_modules.group_to_members[group][members.str_to_int[member]].set();
        result_modules.member_to_groups[member][groups.str_to_int[group]].set();
    }

    file_modules.close();

    return {result_modules, groups, members};
}
