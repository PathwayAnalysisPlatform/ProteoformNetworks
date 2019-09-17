#ifndef PHEGENI_HPP_
#define PHEGENI_HPP_

#include <fstream>
#include <string>
#include <bitset>
#include <unordered_map>
#include <cstring>
#include <iostream>

#include "bimap_str_int.hpp"
#include "conversions.hpp"
#include "networks.hpp"

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

struct load_phegeni_genes_and_traits_result {
    bimap_str_int phegeni_genes;
    bimap_str_int phegeni_traits;
};

load_phegeni_genes_and_traits_result loadPheGenIGenesAndTraits(
        std::string_view path_file_phegeni,
        const bimap_str_int &acceptable_genes);

trait_modules loadPheGenIGeneModules(std::string_view path_file_phegeni,
                                     const bimap_str_int &genes,
                                     const bimap_str_int &traits);

//uss getGeneStrings(const std::bitset<PHEGENI_GENES>& gene_set, const bimap_str_int& genes);
//uss getProteinStrings(const std::bitset<PHEGENI_PROTEINS>& protein_set, const bimap_str_int& proteins);
//uss getProteoformStrings(const std::bitset<PHEGENI_PROTEOFORMS>& proteoform_set, const bimap_str_int& proteoforms);
//
trait_modules convertModulesWithMapping(
        const trait_modules &original_modules,
        const bimap_str_int &original_entities,
        const bimap_str_int &destination_entities,
        const bimap_str_int &traits,
        const ummss &mapping);

trait_modules createPheGenIProteinModules(const trait_modules &gene_modules,
                                          const bimap_str_int &genes,
                                          const bimap_str_int &proteins,
                                          const bimap_str_int &traits,
                                          std::string_view path_file_proteins_to_genes,
                                          std::string_view path_file_protein_edges);

trait_modules removeDisconnectedMembers(trait_modules modules, ummii interactions);

#endif // !PHEGENI_HPP_

