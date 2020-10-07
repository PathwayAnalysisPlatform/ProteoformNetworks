#ifndef PROTEOFORMNETWORKS_NETWORKS_HPP
#define PROTEOFORMNETWORKS_NETWORKS_HPP

#include <string_view>
#include "types.hpp"
#include "bimap_str_int.hpp"
#include <iomanip>
#include <list>
#include "overlap_types.hpp"
#include "maps.hpp"
#include "../config.hpp"

using index_adj_list = vusi;

index_adj_list loadInteractionNetwork(std::string_view path_file_interactions,
                                      const bimap_str_int &entities,
                                      bool has_header_row = true);

struct load_modules_result {
    All_modules entity_modules;
    bimap_str_int groups;
    bimap_str_int members;
};

load_modules_result loadModules(std::string_view path_file_modules,
                                const bimap_str_int &groups, const bimap_str_int &members,
                                bool has_header = true);

load_modules_result loadModules(std::string_view path_file_modules,
                                bool has_header = true);


std::string get_file_name_for_module(std::string module_name);

// Creates a file with all the All_modules: to read all All_modules at once
void writeModulesSingleFile(std::string_view path_file_modules,
                            const All_modules &entity_modules,
                            const bimap_str_int &groups,
                            const bimap_str_int &members
);

// Store All_modules of a level in files, one file for each. For fast access in other python functions.
void writeModulesManyFiles(std::string_view path_file_modules, std::string_view level, std::string_view suffix,
                           const All_modules &entity_modules,
                           const bimap_str_int &groups,
                           const bimap_str_int &members,
                           const vusi &interactions);

All_modules removeDisconnectedMembers(All_modules &modules, const bimap_str_int &groups, const bimap_str_int &members,
                                      const vusi &interactions);

// Create a file with the sizes of the All_modules for each trait at a single level (genes, proteins or proteoforms)
std::map<const std::string, um<int, int>> calculate_and_report_sizes(std::string_view path_reports,
                                                                     const std::map<const std::string, const All_modules> &all_modules,
                                                                     const bimap_str_int &groups);

// Counts the number of edges crossing from one module to the other.
// This function requires that the All_modules vertices and edges files are created
int countCrossingEdges(const std::string &trait1, const std::string &trait2);

#endif //PROTEOFORMNETWORKS_NETWORKS_HPP
