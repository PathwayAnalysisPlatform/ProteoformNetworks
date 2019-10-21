#include <iomanip>
#include <list>
#include "networks.hpp"
#include "overlap_types.hpp"
#include "maps.hpp"

// Read entity interactions in Reactome from a PathwayMatcher edges file
// For each interacting entities A and B, adds both the edges A --> B and B --> A
// The file must have THREE or MORE columns separated by a tab ('\t'). The first two columns contain the source
// and destination interactors. From the third column there are other attributes of the interaction.
vusi loadInteractionNetwork(std::string_view path_file_interactions,
                            const bimap_str_int &entities,
                            bool has_header_row) {

    vusi interactions(entities.int_to_str.size(), std::unordered_set<int>());
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
        if (hasKey(entities.str_to_int, e1) && hasKey(entities.str_to_int, e2)) {
            auto index_e1 = entities.str_to_int.at(e1);
            auto index_e2 = entities.str_to_int.at(e2);
            interactions[index_e1].insert(index_e2);
            interactions[index_e2].insert(index_e1);
        } else {
//            if (!hasKey(entities.str_to_int, e1)) {
//                std::cerr << "Not found entity: **" << e1 << "**" << std::endl;
//            }
//            if (!hasKey(entities.str_to_int, e2)) {
//                std::cerr << "Not found entity: **" << e2 << "**" << std::endl;
//            }
        }
    }

    return interactions;
}

modules removeDisconnectedMembers(modules &modules, const bimap_str_int &groups, const bimap_str_int &members,
                                  const vusi &interactions) {
    // For each module:
    for (auto group_entry = modules.group_to_members.begin();
         group_entry != modules.group_to_members.end(); group_entry++) {
        // Get module members
        int group_index = groups.str_to_int.at(group_entry->first);
        std::unordered_set<int> member_indexes;

        for (int I = 0; I < members.int_to_str.size(); I++)
            if (group_entry->second[I])
                member_indexes.insert(I);

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
                group_entry->second[member_index].reset();

                // Remove owner from vertex owners
                std::string member = members.int_to_str[member_index];
                modules.member_to_groups[member][group_index].reset();
            }
        }

    }

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
        result_modules.group_to_members.emplace(group, base::dynamic_bitset<>(members.int_to_str.size()));
    for (const auto &member : members.int_to_str)
        result_modules.member_to_groups.emplace(member, base::dynamic_bitset<>(groups.int_to_str.size()));

    // Set the members of each group and the owners of each member
    if (has_header) {
        std::getline(file_modules, line);
    }
    while (std::getline(file_modules, group, '\t')) {
        std::getline(file_modules, member);
        result_modules.group_to_members[group][members.str_to_int[member]].set();
        result_modules.member_to_groups[member][groups.str_to_int[group]].set();
    }
    file_modules.close();

    return {result_modules, groups, members};
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

// Creates a file with all the modules: to read all modules at once
void writeModulesSingleFile(std::string_view path_file_modules, std::string_view level, std::string_view suffix,
                            const modules &entity_modules, const bimap_str_int &groups,
                            const bimap_str_int &members) {
    std::string all_traits_file_path = path_file_modules.data();
    all_traits_file_path += level;
    all_traits_file_path += "_modules";
    all_traits_file_path += suffix;

    std::ofstream file_all_traits(all_traits_file_path); // File for all modules
    if (!file_all_traits.is_open()) {
        throw std::runtime_error("Problem opening modules file for writing.\n");
    }
    file_all_traits << "TRAIT\tMEMBER\n";
    for (const auto &group_entry : entity_modules.group_to_members) {

        for (int I = 0; I < members.int_to_str.size(); I++) {
            if (group_entry.second[I]) {
                file_all_traits << group_entry.first << "\t" << members.int_to_str[I] << "\n";
            }
        }
    }
    file_all_traits.close();
}

// Store modules of a level in files, one file for each. For fast access in other python functions.
// For each Trait, it creates two files: vertices (member entities) and edges file (interactions)
void writeModulesManyFiles(std::string_view path_file_modules, std::string_view level, std::string_view suffix,
                           const modules &entity_modules,
                           const bimap_str_int &groups,
                           const bimap_str_int &members,
                           const vusi &interactions) {
    for (const auto &group_entry : entity_modules.group_to_members) {
        std::string file_path_single_trait_vertices = path_file_modules.data();
        std::string file_path_single_trait_edges;
        file_path_single_trait_vertices += get_file_name_for_module(group_entry.first);
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
        for (int I = 0; I < members.int_to_str.size(); I++) {
            if (group_entry.second[I]) {
                file_single_trait_vertices << members.int_to_str[I] << '\n';

                // Write the neighbours of each member.
                // Writes only the edges which go from a lower index to a higher index
                for (int neighbor : interactions[I]) {
                    if (I < neighbor)
                        file_single_trait_edges << members.int_to_str[I] << '\t'
                                                << members.int_to_str[neighbor] << '\n';
                }
            }
        }
        file_single_trait_edges.close();
        file_single_trait_vertices.close();
    }
}

// Calculates and create a file with the sizes of all trait modules at a level (genes, proteins or proteoforms)
um<std::string, int>
calculate_and_report_sizes(std::string_view path_reports, std::string level, modules entity_modules) {
    um<std::string, int> sizes;
    std::ofstream output(path_reports.data() + static_cast<std::string>("module_sizes_") + level + ".tsv");

    if (!output.is_open()) {
        std::string message = "Cannot open report file at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    output << "MODULE\tSIZE\n";
    for (const auto &module_entry : entity_modules.group_to_members) {
        output << module_entry.first << "\t" << module_entry.second.count() << "\n";
        sizes[module_entry.first] = module_entry.second.count();
    }
    output.close();
    return sizes;
}
