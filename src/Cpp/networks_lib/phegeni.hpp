#ifndef PHEGENI_HPP_
#define PHEGENI_HPP_

# include <fstream>
#include <string>
#include <bitset>
#include <unordered_map>
#include <cstring>
#include <iostream>

#include "bimap_str_int.hpp"
#include "conversions.hpp"
#include "networks.hpp"
#include "overlap_types.hpp"
#include "maps.hpp"

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

module_bimaps loadPheGenIGenesAndTraits(
        std::string_view path_file_phegeni,
        const bimap_str_int &acceptable_genes);

// Read Phegeni trait modules with genes as members. Only gene members also in the acceptable gene list.
// The only method to create gene modules. For proteins and proteoforms use the method: createPheGenIModules
modules loadPheGenIGeneModules(std::string_view path_file_phegeni,
                               const bimap_str_int &genes,
                               const bimap_str_int &traits);

// Creates new modules with the same groups, but replacing the member identifiers to new identifiers with the mapping
modules convertModulesWithMapping(
        const modules &original_modules,
        const bimap_str_int &original_entities,
        const bimap_str_int &destination_entities,
        const bimap_str_int &traits,
        const ummss &mapping);

// Creates modules of proteins and proteoforms.
// Converts the modules from gene --> protein or protein --> proteoform, and then removes the disconnected vertices
modules createPheGenIModules(const modules &prev_modules,
                             const bimap_str_int &prev_entities,
                             const bimap_str_int &entities,
                             const bimap_str_int &traits,
                             std::string_view path_file_mapping_prev_entities_to_entities,
                             std::string_view path_file_entity_interactions);

#endif // !PHEGENI_HPP_

