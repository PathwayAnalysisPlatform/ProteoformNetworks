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


// Create a file with the sizes of the All_modules for each trait at a single level (genes, proteins or proteoforms)
std::map<const std::string, um<int, int>> calculate_and_report_sizes(std::string_view path_reports,
                                                                     const std::map<const std::string, const All_modules> &all_modules,
                                                                     const bimap_str_int &groups);

// Counts the number of edges crossing from one module to the other.
// This function requires that the All_modules vertices and edges files are created
int countCrossingEdges(const std::string &trait1, const std::string &trait2);

#endif //PROTEOFORMNETWORKS_NETWORKS_HPP
