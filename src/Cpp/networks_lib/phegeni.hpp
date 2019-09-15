#ifndef PHEGENI_HPP

#include <fstream>
#include <string>
#include <bitset>
#include <unordered_map>

#include "bimap_str_int.hpp"

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

const size_t PHEGENI_TRAITS = 790;
const size_t PHEGENI_GENES = 3350;
const size_t PHEGENI_PROTEINS = 5695;
const size_t PHEGENI_PROTEOFORMS = 6971;

struct load_phegeni_genes_and_traits_result {
    bimap_str_int phegeni_genes;
    bimap_str_int phegeni_traits;
};

struct trait_modules {
    msb traits_to_genes;
    msb genes_to_traits;
};

load_phegeni_genes_and_traits_result loadPheGenIGenesAndTraits(
        std::string_view path_file_phegeni,
        const bimap_str_int &acceptable_genes);

trait_modules loadPheGenIModules(
        std::string_view path_file_phegeni,
        const bimap_str_int &acceptable_genes);

umsb convertGeneSets(
        const umsb &traits_to_genes,
        const bimap_str_int &phegeni_genes,
        const ummss &mapping_genes_to_proteins,
        const bimap_str_int &proteins,
        const ummss &adjacency_list_proteins);

umsb convertProteinSets(const umsb &traits_to_proteins,
                        const bimap_str_int &proteins,
                        const ummss &mapping_proteins_to_proteoforms,
                        const bimap_str_int &proteoforms,
                        const ummss &adjacency_list_proteoforms);

//uss getGeneStrings(const std::bitset<PHEGENI_GENES>& gene_set, const bimap_str_int& genes);
//uss getProteinStrings(const std::bitset<PHEGENI_PROTEINS>& protein_set, const bimap_str_int& proteins);
//uss getProteoformStrings(const std::bitset<PHEGENI_PROTEOFORMS>& proteoform_set, const bimap_str_int& proteoforms);
//
umsb convertSets(
        const umsb &traits_to_original_entities,
        const vs &index_to_original_entities,
        const ummss &mapping,
        const umsi &result_entities_to_index,
        const ummss &adjacency_list_result_entities);


#endif // !PHEGENI_HPP

