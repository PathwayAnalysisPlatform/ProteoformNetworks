#ifndef PHEGENI_HPP
#define PHEGENI_HPP

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

using phegeni_gene_to_traits = um<std::string, std::bitset<PHEGENI_TRAITS>>;	// Gene -> Traits where it participates
using phegeni_trait_to_genes = um<std::string, std::bitset<PHEGENI_GENES>>;
using phegeni_trait_to_proteins = um<std::string, std::bitset<PHEGENI_PROTEINS>>;
using phegeni_trait_to_proteoforms = um<std::string, std::bitset<PHEGENI_PROTEOFORMS>>;

struct load_phegeni_genes_and_traits_result {
	const bimap_str_int phegeni_genes;
	const bimap_str_int phegeni_traits;
};

load_phegeni_genes_and_traits_result loadPheGenIGenesAndTraits(
	std::string_view path_file_phegeni,
	const bimap_str_int& reactome_genes);

struct load_phegeni_sets_result {
	phegeni_trait_to_genes traits_to_genes;
	phegeni_gene_to_traits genes_to_traits;
};

load_phegeni_sets_result loadPheGenISets(
	std::string_view path_file_phegeni,
	const bimap_str_int& reactome_genes,
	const bimap_str_int& phegeni_genes,
	const bimap_str_int& phegeni_traits);

phegeni_trait_to_proteins convertGeneSets(
	const phegeni_trait_to_genes& traits_to_genes,
	const bimap_str_int& phegeni_genes,
	const ummss& mapping_genes_to_proteins,
	const bimap_str_int& proteins,
	const ummss& adjacency_list_proteins);

std::unordered_map<std::string, std::bitset<PHEGENI_PROTEOFORMS>> convertProteinSets(const std::unordered_map<std::string, std::bitset<PHEGENI_PROTEINS>>& traits_to_proteins,
	const bimap_str_int& proteins,
	const ummss& mapping_proteins_to_proteoforms,
	const bimap_str_int& proteoforms,
	const ummss& adjacency_list_proteoforms);

//uss getGeneStrings(const std::bitset<PHEGENI_GENES>& gene_set, const bimap_str_int& genes);
//uss getProteinStrings(const std::bitset<PHEGENI_PROTEINS>& protein_set, const bimap_str_int& proteins);
//uss getProteoformStrings(const std::bitset<PHEGENI_PROTEOFORMS>& proteoform_set, const bimap_str_int& proteoforms);
//
umss createTraitNames(const um<std::string, std::bitset<PHEGENI_GENES>>& traits_to_genes);

template <size_t total_num_original_entities, size_t total_num_result_entities>
um<std::string, std::bitset<total_num_result_entities>> convertSets(
	const um<std::string, std::bitset<total_num_original_entities>>& traits_to_original_entities,
	const vs& index_to_original_entities,
	const ummss& mapping,
	const umsi& result_entities_to_index,
	const ummss& adjacency_list_result_entities) {
	std::unordered_map<std::string, std::bitset<total_num_result_entities>> traits_to_result_entities;
	for (const auto& trait_entry : traits_to_original_entities) {  // For each trait entry
		uss candidates;
		for (int I = 0; I < total_num_original_entities; I++) {  // For each member in this trait
			if (trait_entry.second.test(I)) {
				auto range = mapping.equal_range(index_to_original_entities[I]);  // For each mapping entity
				for (auto it = range.first; it != range.second; it++) {
					candidates.insert(it->second);
				}
			}
		}

		// Keep only those connected to any of the other gene set members in the reference network
		for (const auto& candidate : candidates) {
			auto range = adjacency_list_result_entities.equal_range(candidate);
			for (auto it = range.first; it != range.second; ++it) {
				if (candidates.find(it->second) != candidates.end()) {
					traits_to_result_entities[trait_entry.first].set(result_entities_to_index.at(candidate));  // Set in stone that the candidate is in the new set
					break;
				}
			}
		}
	}

	return traits_to_result_entities;
}

#endif // !PHEGENI_HPP

