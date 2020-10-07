#include "networks.hpp"
#include "../config.hpp"

// Read entity interactions file.
// For each interacting entities A and B, adds both the edges A --> B
// The file is a tab separated file, with two columns: Interactor A and Interactor B. Anything after the second column
// is discarded.
index_adj_list loadInteractionNetwork(std::string_view path_file_interactions,
                                      const bimap_str_int &entities,
                                      bool has_header_row) {

    index_adj_list interactions(entities.size(), std::unordered_set<int>());
    std::ifstream file_interactions(path_file_interactions.data());
    std::string A, B, other_fields;

    if (!file_interactions.is_open()) {
        std::string message = "Error reading file at: ";
        message += __FUNCTION__;
        throw std::runtime_error(message);
    }

    if (has_header_row)
        getline(file_interactions, other_fields);    // Discard the header line.
    while (getline(file_interactions, A, '\t')) {
        getline(file_interactions, B, '\t');   // get second entity
        getline(file_interactions, other_fields);    // discard remaining characters of the line
        if (entities.has(A) && entities.has(B)) {
            auto index_A = entities.index_of(A);
            auto index_B = entities.index_of(B);
            interactions[index_A].insert(index_B);
            interactions[index_B].insert(index_A);
        } else {
            if (!entities.has(A)) {
                std::cerr << "Not found entity: **" << A << "**" << std::endl;
            }
            if (!entities.has(B)) {
                std::cerr << "Not found entity: **" << B << "**" << std::endl;
            }
        }
    }

    return interactions;
}

All_modules removeDisconnectedMembers(All_modules &modules,
                                      const bimap_str_int &groups, const bimap_str_int &members,
                                      const vusi &interactions) {
    // For each module:
    for (int group_index = 0; group_index < modules.group_to_members.size(); group_index++) {
        // Get module members
        std::unordered_set<int> member_indexes;

        for (int member_index = 0; member_index < members.int_to_str.size(); member_index++)
            if (modules.group_to_members[group_index][member_index])
                member_indexes.insert(member_index);

        // For each member, iterate over its interactors, check if any is also member of the group
        for (int member_index : member_indexes) {
            bool isConnected = false;

            // Check if any interactor of this member is in the group
            for (const auto &interactor : interactions[member_index]) {
                if (hasValue(member_indexes, interactor)) {   // If the interactor is in the group
                    isConnected = true;
                    break;
                }
            }

            if (!isConnected) {
                // Remove vertex from the group
                modules.group_to_members[group_index][member_index].reset();

                // Remove owner from vertex owners
                modules.member_to_groups[member_index][group_index].reset();
            }
        }

    }

    return modules;
}

// Reads groups from files.
// The file should follow the format convention: two colums, no extra spaces at the end of each row.
// Columns separated by a '\t' (tab) character.
load_modules_result loadModules(std::string_view path_file_modules,
                                const bimap_str_int &groups, const bimap_str_int &members,
                                bool has_header) {
    All_modules result_modules;
    std::string line, group, member;

    // # Read the file to create the All_modules with its members
    std::ifstream file_modules(path_file_modules.data());
    if (!file_modules.is_open()) {
        std::string message = "Cannot open path_file_modules at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    // ## Read All_modules
    // Initialize the All_modules to the correct sizes
    for (int group_index = 0; group_index < groups.size(); group_index++)
        result_modules.group_to_members.push_back(base::dynamic_bitset<>(members.size()));
    for (int member_index = 0; member_index < members.size(); member_index++)
        result_modules.member_to_groups.push_back(base::dynamic_bitset<>(groups.size()));

    // Set the members of each group and the owners of each member
    if (has_header) {
        std::getline(file_modules, line);
    }
    while (std::getline(file_modules, group, '\t')) {
        std::getline(file_modules, member);
        try {
            int group_index = groups.str_to_int.at(group);
            int member_index = members.str_to_int.at(member);
            result_modules.group_to_members[group_index][member_index].set();
            result_modules.member_to_groups[member_index][group_index].set();
        }
        catch (const std::out_of_range &oor) {
            std::cerr << "Out of Range error: " << oor.what() << '\n';
        }
    }
    file_modules.close();

    return {result_modules, groups, members};
}

// Reads groups from files. It creates the bimap structures for the groups and the members from the file.
// The file should follow the format convention: two colums, no extra spaces at the end of each row.
// Columns separated by a '\t' (tab) character.
load_modules_result loadModules(std::string_view path_file_modules, bool has_header) {
    All_modules result_modules;
    bimap_str_int groups, members;

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

    return loadModules(path_file_modules, groups, members, has_header);
}

std::string get_file_name_for_module(std::string module_name) {
    // Remove all quote ('"') occurrences
    std::list<char> unwanted_chars = {'"', ' ', ','};
    for (char c : unwanted_chars) {
        size_t position = 0;
        while (true) {
            position = module_name.find(c, position);
            if (position == std::string::npos) {
                break;
            }
            module_name.replace(position, 1, c == '"' ? "" : "_");
        }
    }
    return module_name;
}

// Creates a file with all the All_modules: to read all All_modules at once
// The format of the file is a table with two columns. The lines are records of tuples trait-member
void writeModulesSingleFile(std::string_view path_file_modules,
                            const All_modules &entity_modules,
                            const bimap_str_int &groups,
                            const bimap_str_int &members) {
    std::cerr << "Writing All_modules to file: " << path_file_modules << std::endl;
    std::ofstream file_all_traits(path_file_modules.data()); // File for all All_modules
    if (!file_all_traits.is_open()) {
        throw std::runtime_error("Problem opening All_modules file for writing.\n");
    }
    file_all_traits << "TRAIT\tMEMBER\n";
    for (int g = 0; g < groups.size(); g++) {
        for (int m = 0; m < members.size(); m++) {
            if (entity_modules.group_to_members[g][m]) {
                file_all_traits << groups.int_to_str[g] << "\t" << members.int_to_str[m] << "\n";
            }
        }
    }
    file_all_traits.close();
}

// Store All_modules of a level, each module in separate files. For fast access in other python functions.
// For each Trait, it creates two files: vertices (member entities) and edges file (interactions)
void writeModulesManyFiles(std::string_view path_file_modules,
                           std::string_view level,
                           std::string_view suffix,
                           const All_modules &entity_modules,
                           const bimap_str_int &groups,
                           const bimap_str_int &members,
                           const vusi &interactions) {
    for (int g = 0; g < groups.size(); g++) {
        std::string file_path_single_trait_vertices = path_file_modules.data();
        std::string file_path_single_trait_edges;
        file_path_single_trait_vertices += get_file_name_for_module(groups.int_to_str[g]);
        file_path_single_trait_vertices += "_";
        file_path_single_trait_vertices += level.data();
        file_path_single_trait_vertices += "_";

        file_path_single_trait_edges = file_path_single_trait_vertices;

        file_path_single_trait_vertices += "vertices";
        file_path_single_trait_vertices += suffix;

        file_path_single_trait_edges += "edges";
        file_path_single_trait_edges += suffix;

        std::ofstream file_single_trait_vertices(file_path_single_trait_vertices);
        std::ofstream file_single_trait_edges(file_path_single_trait_edges);

        file_single_trait_vertices << "MEMBER\n";
        file_single_trait_edges << "MEMBER1\tMEMBER2\n";
        for (int m = 0; m < members.size(); m++) {
            if (entity_modules.group_to_members[g][m]) {
                file_single_trait_vertices << members.int_to_str[m] << '\n';

                // Write the neighbours of each member.
                // Writes only the edges which go from a lower index to a higher index
                for (int neighbor : interactions[m]) {
                    // If the neighbor is also in the module
                    if (entity_modules.group_to_members[g][neighbor] && m < neighbor)
                        file_single_trait_edges << members.int_to_str[m] << '\t'
                                                << members.int_to_str[neighbor] << '\n';
                }
            }
        }
        file_single_trait_edges.close();
        file_single_trait_vertices.close();
    }
}

// Calculates and create a file with the sizes of all trait All_modules at all levels (genes, proteins or proteoforms)
std::map<const std::string, um<int, int>> calculate_and_report_sizes(std::string_view path_reports,
                                                                     const std::map<const char*, const All_modules> &all_modules,
                                                                     const bimap_str_int &groups) {
    std::map<const std::string, um<int, int>> sizes;

    for (const auto &level : config::LEVELS) {
        std::string file_name = path_reports.data() + static_cast<std::string>("module_sizes_") + level + ".tsv";
        std::ofstream output(file_name);

        if (!output.is_open()) {
            std::string message = "Cannot open report file at ";
            std::string function = __FUNCTION__;
            throw std::runtime_error(message + function);
        }

        std::cerr << "Writing file: " << file_name << std::endl;

        output << "MODULE\tSIZE\n";
        for (int group_index = 0; group_index < all_modules.at(level).group_to_members.size(); group_index++) {
            output << groups.int_to_str[group_index] << "\t"
                   << all_modules.at(level).group_to_members[group_index].count()
                   << "\n";
            sizes[level][group_index] = all_modules.at(level).group_to_members[group_index].count();
        }
        output.close();
    }

    return sizes;
}
