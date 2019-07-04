#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <algorithm>
#include <bitset>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "bimap.hpp"
#include "types.hpp"
#include "entity.hpp"
#include "utility.hpp"
#include "phegeni.hpp"

const size_t NUM_GENES = 23970;
const size_t NUM_PROTEINS = 10778;
const size_t NUM_PROTEOFORMS = 13911;

const size_t NUM_PHEGEN_TRAITS = 846;
const size_t MAX_PHEGEN_TRAITS = 1174;
const size_t NUM_PHEGEN_GENES = 8153;
const size_t MAX_PHEGEN_GENES = 17128;
const size_t NUM_PHEGEN_PROTEINS = 5695;
const size_t NUM_PHEGEN_PROTEOFORMS = 6971;

const int MIN_OVERLAP_SIZE = 1;
const int MAX_OVERLAP_SIZE = 100;

const int MIN_SET_SIZE = 1;
const int MAX_SET_SIZE = 200;

void doOverlapAnalysis(std::string_view path_file_PheGenI, std::string_view path_file_gene_search);

struct Frequencies {
	msi modifications;
	msi proteins;
	msi proteoforms;
};

struct load_trait_gene_sets_result {
	um<std::string, std::bitset<NUM_PHEGEN_TRAITS>> genes_to_sets;
	um<std::string, std::bitset<NUM_PHEGEN_GENES>> sets_to_genes;
};

struct load_mapping_result {
	ummss ones_to_others;
	ummss others_to_ones;
};

umss loadPathwayNames(std::string_view path_protein_search_file);

msi fillMap(const vs& index_to_entities);

std::unordered_map<std::string, std::bitset<NUM_GENES>> loadGeneSets(std::string_view file_path,
	const umsi& entities_to_index,
	bool pathways);
std::unordered_map<std::string, std::bitset<NUM_PROTEINS>> loadProteinSets(std::string_view file_path,
	const umsi& entities_to_index,
	bool pathways);
std::unordered_map<std::string, std::bitset<NUM_PROTEOFORMS>> loadProteoformSets(std::string_view file_path,
	const umsi& entities_to_index,
	bool pathways);

uss getGeneStrings(const std::bitset<NUM_GENES>& gene_set,
	const vs& index_to_genes);
uss getProteinStrings(const std::bitset<NUM_PROTEINS>& protein_set,
	const vs& index_to_proteins);
uss getProteoformStrings(const std::bitset<NUM_PROTEOFORMS>& proteoform_set,
	const vs& index_to_proteoforms);

uss getGeneStrings(const std::bitset<NUM_PHEGEN_GENES>& gene_set,
	const vs& index_to_genes);
uss getProteinStrings(const std::bitset<NUM_PHEGEN_PROTEINS>& protein_set,
	const vs& index_to_proteins);
uss getProteoformStrings(const std::bitset<NUM_PHEGEN_PROTEOFORMS>& proteoform_set,
	const vs& index_to_proteoforms);

load_phegeni_entities_result loadPheGenIEntities(std::string_view path_file_PheGenI,
	const double& max_p_value,
	const umsi& reactome_entities_to_index);

load_trait_gene_sets_result loadTraitGeneSets(const std::string& path_file_phegen,
	const double& max_p_value,
	const bimap& genes,
	const bimap& traits,
	const umsi& reactome_genes_to_index);

load_mapping_result loadMapping(const std::string& path_file_mapping);

std::unordered_map<std::string, std::bitset<NUM_PHEGEN_PROTEINS>> convertGeneSets(const std::unordered_map<std::string, std::bitset<NUM_PHEGEN_GENES>>& traits_to_genes,
	const vs& index_to_genes,
	const ummss& mapping_genes_to_proteins,
	const umsi& proteins_to_index,
	const ummss& adjacency_list_proteins);

std::unordered_map<std::string, std::bitset<NUM_PHEGEN_PROTEOFORMS>> convertProteinSets(const std::unordered_map<std::string, std::bitset<NUM_PHEGEN_PROTEINS>>& traits_to_proteins,
	const vs& index_to_proteins,
	const ummss& mapping_proteins_to_proteoforms,
	const umsi& proteoforms_to_index,
	const ummss& adjacency_list_proteoforms);

template <size_t total_num_entities>
void printMembers(std::ostream& output, const std::bitset<total_num_entities>& entity_set, const vs& index_to_entities) {
	int printed = 0;
	int total = entity_set.count();
	output << "[";
	for (int I = 0; I < total_num_entities; I++) {
		if (entity_set.test(I)) {
			output << "\"" << index_to_entities[I] << "\"";
			printed++;
			if (printed != total) {
				output << ",";
			}
		}
	}
	output << "]";
}

// Version with set size and overlap size limits
template <size_t total_num_entities>
std::map<std::pair<std::string, std::string>, std::bitset<total_num_entities>> findOverlappingPairs(const std::unordered_map<std::string, std::bitset<total_num_entities>>& sets_to_members,
	const int& min_overlap, const int& max_overlap,
	const int& min_set_size, const int& max_set_size) {
	std::map<std::pair<std::string, std::string>, std::bitset<total_num_entities>> result;
	std::vector<typename std::unordered_map<std::string, std::bitset<total_num_entities>>::const_iterator> nav;

	for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
		int set_size = it->second.count();
		if (min_set_size <= set_size && set_size <= max_set_size) {
			nav.push_back(it);
		}
	}

	std::cerr << "elementos: " << nav.size() << ", parejas: " << (nav.size() * nav.size() - nav.size()) / 2 << "\n";
	auto t0 = clock();

	for (auto vit1 = nav.begin(); vit1 != nav.end(); vit1++) {
		for (auto vit2 = vit1 + 1; vit2 != nav.end(); vit2++) {
			std::bitset<total_num_entities> overlap = (*vit1)->second & (*vit2)->second;
			int overlap_size = overlap.count();
			if (min_overlap <= overlap_size && overlap_size <= max_overlap) {
				result.emplace(make_pair((*vit1)->first, (*vit2)->first), overlap);
			}
		}
	}
	auto t1 = clock();
	std::cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
	return result;
}

// This version includes modified ratio of the proteoform members of the set, and another ratio for the overlapping proteoforms.
template <size_t S>
std::map<std::pair<std::string, std::string>, std::bitset<S>> findOverlappingProteoformSets(const std::unordered_map<std::string, std::bitset<S>>& sets_to_members,
	const int& min_overlap, const int& max_overlap,
	const int& min_set_size, const int& max_set_size,
	const std::bitset<S>& modified_proteoforms,
	const double& min_all_modified_ratio,
	const double& min_overlap_modified_ratio) {
	std::vector<typename std::unordered_map<std::string, std::bitset<S>>::const_iterator> nav;
	for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
		int set_size = it->second.count();
		if (min_set_size <= set_size && set_size <= max_set_size) {
			//  cerr << "MIN_SET_SIZE: " << min_set_size << " SET_SIZE: " << set_size << " MAX_SET_SIZE: " << max_set_size << "\n";
			float percentage = static_cast<float>((modified_proteoforms & it->second).count()) / static_cast<float>(set_size);
			if (percentage >= min_all_modified_ratio) {
				nav.push_back(it);
			}
		}
	}

	std::map<std::pair<std::string, std::string>, std::bitset<S>> result;
	std::cerr << "elementos: " << nav.size() << ", parejas: " << (nav.size() * nav.size() - nav.size()) / 2 << "\n";
	auto t0 = clock();
	for (auto vit1 = nav.begin(); vit1 != nav.end(); vit1++) {
		for (auto vit2 = next(vit1); vit2 != nav.end(); vit2++) {
			std::bitset<S> overlap = (*vit1)->second & (*vit2)->second;
			if (min_overlap <= overlap.count() && overlap.count() <= max_overlap) {
				float percentage = static_cast<float>((modified_proteoforms & overlap).count()) / static_cast<float>(overlap.count());
				if (percentage >= min_overlap_modified_ratio) {
					result.emplace(make_pair((*vit1)->first, (*vit2)->first), overlap);
				}
			}
		}
	}
	auto t1 = clock();
	std::cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";

	return result;
}

umss createTraitNames(const um<std::string, std::bitset<NUM_PHEGEN_GENES>>& traits_to_genes);

struct reactome_networks {
	ummss adjacency_list_proteins;
	ummss adjacency_list_proteoforms;
};

reactome_networks loadReactomeNetworks(std::string_view path_file_gene_search,
	std::string_view path_file_protein_search,
	std::string_view path_file_proteoform_search);

ummss loadProteinsAdjacencyList(std::string_view search_file_path);

ummss loadProteoformsAdjacencyList(std::string_view search_file_path);

bimap deductProteinsFromGenes(std::string_view path_file_mapping_proteins_genes,
	const umsi& genes_to_index,
	const ummss& genes_to_proteins);

bimap deductProteoformsFromProteins(const ummss& proteins_to_proteoforms, const umsi& proteins_to_index);



#endif /* OVERLAP_H_ */