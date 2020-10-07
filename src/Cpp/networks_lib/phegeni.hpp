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
#include "../config.hpp"
#include "../overlap_analysis.hpp"

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

bimaps loadPheGenIGenesAndTraits(
        std::string_view path_file_phegeni,
        const bimap_str_int &acceptable_genes);

All_modules get_or_create_gene_modules(std::string_view path_file_phegeni,
                                       const bimap_str_int &traits,
                                       const entities_bimap &genes,
                                       const entities_bimap &small_molecules,
                                       std::string_view path_file_gene_interactions,
                                       const std::string &path_output,
                                       bool keep_disconnected_nodes);

// Read Phegeni trait All_modules with genes as members. Only gene members also in the acceptable gene list.
// The only method to create gene All_modules. For proteins and proteoforms use the method: createPheGenIModules
All_modules createGeneModules(std::string_view path_file_phegeni,
                              const bimap_str_int &genes, const bimap_str_int &small_molecules,
                              const bimap_str_int &traits, const std::string &path_file_gene_interactions,
                              const std::string &path_modules, bool keep_disconnected_nodes);

// Creates new All_modules with the same groups, but replacing the member identifiers to new identifiers with the mapping
All_modules convertModulesWithMapping(
        const All_modules &original_modules,
        const bimap_str_int &original_entities,
        const bimap_str_int &destination_entities,
        const bimap_str_int &traits,
        const ummss &mapping);

// Creates All_modules of proteins and proteoforms.
// Converts the All_modules from gene --> protein or protein --> proteoform, and then removes the disconnected vertices
All_modules createProteinOrProteoformModules(const All_modules &prev_modules, const bimap_str_int &prev_entities,
                                             const bimap_str_int &entities, const bimap_str_int &traits,
                                             const ummss &mapping,
                                             const std::string &path_file_entity_interactions,
                                             const std::string &path_modules,
                                             const std::string &level, const std::string &suffix,
                                             bool keep_disconnected_nodes);

#endif // !PHEGENI_HPP_

