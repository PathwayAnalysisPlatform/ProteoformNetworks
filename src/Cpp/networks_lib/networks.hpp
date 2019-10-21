#ifndef PROTEOFORMNETWORKS_NETWORKS_HPP
#define PROTEOFORMNETWORKS_NETWORKS_HPP

#include <string_view>
#include "types.hpp"
#include "bimap_str_int.hpp"

vusi loadInteractionNetwork(std::string_view path_file_interactions,
                            const bimap_str_int &entities,
                            bool has_header_row = true);

struct load_modules_result {
    modules entity_modules;
    bimap_str_int groups;
    bimap_str_int members;
};

load_modules_result loadModules(std::string_view path_file_modules, bool has_header = true);

std::string get_file_name_for_module(std::string module_name);

// Creates a file with all the modules: to read all modules at once
void writeModulesSingleFile(std::string_view path_file_modules, std::string_view level, std::string_view suffix,
                            const modules &entity_modules,
                            const bimap_str_int &groups,
                            const bimap_str_int &members
);

// Store modules of a level in files, one file for each. For fast access in other python functions.
void writeModulesManyFiles(std::string_view path_file_modules, std::string_view level, std::string_view suffix,
                           const modules &entity_modules,
                           const bimap_str_int &groups,
                           const bimap_str_int &members,
                           const vusi &interactions);

modules removeDisconnectedMembers(modules &modules, const bimap_str_int &groups, const bimap_str_int &members,
                                  const vusi &interactions);

// Create a file with the sizes of the modules for each trait at a single level (genes, proteins or proteoforms)
um<std::string, int> calculate_and_report_sizes(std::string_view path_reports,
                                                std::string level, modules entity_modules);

#endif //PROTEOFORMNETWORKS_NETWORKS_HPP
