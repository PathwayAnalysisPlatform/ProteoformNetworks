#ifndef PHEGENI_HPP
#define PHEGENI_HPP

#include <fstream>
#include <string>
#include <bitset>
#include <unordered_map>

#include "bimap.hpp"

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

const size_t PHEGENI_TRAITS = 790;
const size_t PHEGENI_GENES = 3350;
const size_t PHEGENI_PROTEINS = 5695;
const size_t PHEGENI_PROTEOFORMS = 6971;

using  phegeni_gene_to_traits = um<std::string, std::bitset<PHEGENI_TRAITS>>;	// Gene -> Traits where it participates
using phegeni_trait_to_genes = um<std::string, std::bitset<PHEGENI_GENES>>;
using phegeni_protein_to_traits = um<std::string, std::bitset<PHEGENI_PROTEINS>>;
using phegeni_proteoform_to_traits = um<std::string, std::bitset<PHEGENI_PROTEOFORMS>>;

struct load_phegeni_entities_result {
	const bimap genes;
	const bimap traits;
};

load_phegeni_entities_result loadPheGenIEntities(
	std::string_view path_file_phegeni,
	const bimap& reactome_genes,
	const double& max_p_value);

struct load_phegeni_sets_result {
	phegeni_trait_to_genes trait_to_genes;
	phegeni_gene_to_traits gene_to_traits;
};

load_phegeni_sets_result loadPheGenISets(
	std::string_view path_file_phegeni,
	const bimap& reactome_genes,
	const double& max_p_value,
	const bimap& phegeni_genes,
	const bimap& phegeni_traits);

//phegeni_protein_sets convertGeneSets(
//	const phegeni_gene_sets& traits_to_genes,
//	const bimap& genes,
//	const ummss& mapping_genes_to_proteins,
//	const bimap& proteins,
//	const ummss& adjacency_list_proteins);
//
//std::unordered_map<std::string, std::bitset<PHEGENI_PROTEOFORMS>> convertProteinSets(const std::unordered_map<std::string, std::bitset<PHEGENI_PROTEINS>>& traits_to_proteins,
//	const bimap& proteins,
//	const ummss& mapping_proteins_to_proteoforms,
//	const bimap& proteoforms,
//	const ummss& adjacency_list_proteoforms);

//uss getGeneStrings(const std::bitset<PHEGENI_GENES>& gene_set, const bimap& genes);
//uss getProteinStrings(const std::bitset<PHEGENI_PROTEINS>& protein_set, const bimap& proteins);
//uss getProteoformStrings(const std::bitset<PHEGENI_PROTEOFORMS>& proteoform_set, const bimap& proteoforms);
//
//umss createTraitNames(const um<std::string, std::bitset<PHEGENI_GENES>>& traits_to_genes);

#endif // !PHEGENI_HPP

