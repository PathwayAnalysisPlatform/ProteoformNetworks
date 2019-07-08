#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <iostream>
#include <ctime>

#include "phegeni.hpp"
#include "bimap.hpp"
#include "uniprot.hpp"
#include "reactome.hpp"

const int MIN_OVERLAP_SIZE = 1;
const int MAX_OVERLAP_SIZE = 100;

const int MIN_SET_SIZE = 1;
const int MAX_SET_SIZE = 200;

struct Frequencies {
	msi modifications;
	msi proteins;
	msi proteoforms;
};

void doOverlapAnalysis(
	std::string_view path_file_PheGenI,
	std::string_view path_file_gene_search,
	std::string_view path_file_mapping_proteins_to_genes,
	std::string_view path_file_proteoform_search);

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

template <size_t total_num_entities>
void printMembers(std::ostream& output, const std::bitset<total_num_entities>& entity_set, const bimap& entities) {
	int printed = 0;
	int total = entity_set.count();
	output << "[";
	for (int I = 0; I < total_num_entities; I++) {
		if (entity_set.test(I)) {
			output << "\"" << entities.entities[I] << "\"";
			printed++;
			if (printed != total) {
				output << ",";
			}
		}
	}
	output << "]";
}

void printMembers(std::ostream& output, const uss& members);

template<size_t total_num_entities>
uss getStringSetFromEntityBitset(const std::bitset<total_num_entities>& entity_set, const bimap& entities) {
	uss result;
	for (int I = 0; I < total_num_entities; I++) {
		if (entity_set.test(I)) {
			result.insert(entities.entities[I]);
		}
	}
	return result;
}

#endif /* OVERLAP_H_ */