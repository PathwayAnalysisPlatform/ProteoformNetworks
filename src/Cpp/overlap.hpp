#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <algorithm>
#include <bitset>
#include <charconv>
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

#include "entity.hpp"

const size_t NUM_GENES = 23970;
const size_t NUM_PROTEINS = 10778;
const size_t NUM_PROTEOFORMS = 13911;

const size_t MAX_NUM_PHEGEN_TRAITS = 1174;
const size_t MAX_NUM_PHEGEN_GENES = 17128;
const size_t NUM_PHEGEN_TRAITS = 846;
const size_t NUM_PHEGEN_GENES = 8153;
const size_t NUM_PHEGEN_PROTEINS = 5695;
const size_t NUM_PHEGEN_PROTEOFORMS = 6971;

const int MIN_OVERLAP_SIZE = 1;
const int MAX_OVERLAP_SIZE = 100;

const int MIN_SET_SIZE = 1;
const int MAX_SET_SIZE = 200;

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

struct Frequencies {
   std::map<std::string, int> modifications;
   std::map<std::string, int> proteins;
   std::map<std::string, int> proteoforms;
};

struct index_to_entitites_phegen_result {
   std::vector<std::string> index_to_genes;
   std::vector<std::string> index_to_traits;
};

struct load_genes_phegen_result {
   std::vector<std::string> index_to_genes;
   std::vector<std::string> index_to_traits;
   std::unordered_map<std::string, int> genes_to_index;
   std::unordered_map<std::string, int> traits_to_index;
};

struct load_trait_gene_sets_result {
   std::unordered_map<std::string, std::bitset<NUM_PHEGEN_TRAITS>> genes_to_sets;
   std::unordered_map<std::string, std::bitset<NUM_PHEGEN_GENES>> sets_to_genes;
};

struct load_mapping_result {
   std::unordered_multimap<std::string, std::string> ones_to_others;
   std::unordered_multimap<std::string, std::string> others_to_ones;
};

std::unordered_map<std::string, std::string> loadPathwayNames(std::string_view path_protein_search_file);

std::map<std::string, int> fillMap(const std::vector<std::string>& index_to_entities);

std::unordered_map<std::string, std::bitset<NUM_GENES>> loadGeneSets(std::string_view file_path,
                                                                     const std::unordered_map<std::string, int>& entities_to_index,
                                                                     bool pathways);
std::unordered_map<std::string, std::bitset<NUM_PROTEINS>> loadProteinSets(std::string_view file_path,
                                                                           const std::unordered_map<std::string, int>& entities_to_index,
                                                                           bool pathways);
std::unordered_map<std::string, std::bitset<NUM_PROTEOFORMS>> loadProteoformSets(std::string_view file_path,
                                                                                 const std::unordered_map<std::string, int>& entities_to_index,
                                                                                 bool pathways);

std::unordered_set<std::string> getGeneStrings(const std::bitset<NUM_GENES>& gene_set,
                                               const std::vector<std::string>& index_to_genes);
std::unordered_set<std::string> getProteinStrings(const std::bitset<NUM_PROTEINS>& protein_set,
                                                  const std::vector<std::string>& index_to_proteins);
std::unordered_set<std::string> getProteoformStrings(const std::bitset<NUM_PROTEOFORMS>& proteoform_set,
                                                     const std::vector<std::string>& index_to_proteoforms);

std::unordered_set<std::string> getGeneStrings(const std::bitset<NUM_PHEGEN_GENES>& gene_set,
                                               const std::vector<std::string>& index_to_genes);
std::unordered_set<std::string> getProteinStrings(const std::bitset<NUM_PHEGEN_PROTEINS>& protein_set,
                                                  const std::vector<std::string>& index_to_proteins);
std::unordered_set<std::string> getProteoformStrings(const std::bitset<NUM_PHEGEN_PROTEOFORMS>& proteoform_set,
                                                     const std::vector<std::string>& index_to_proteoforms);

load_genes_phegen_result loadGenesPheGen(std::string_view path_file_PheGenI_full,
                                         const double& max_p_value,
                                         const std::unordered_map<std::string, int>& reactome_entities_to_index);

load_trait_gene_sets_result loadTraitGeneSets(const std::string& path_file_phegen,
                                              const double& max_p_value,
                                              const std::vector<std::string>& index_to_genes,
                                              const std::vector<std::string>& index_to_traits,
                                              const std::unordered_map<std::string, int>& genes_to_index,
                                              const std::unordered_map<std::string, int>& traits_to_index,
                                              const std::unordered_map<std::string, int>& reactome_genes_to_index);

load_mapping_result loadMapping(const std::string& path_file_mapping);

std::unordered_map<std::string, std::bitset<NUM_PHEGEN_PROTEINS>> convertGeneSets(const std::unordered_map<std::string, std::bitset<NUM_PHEGEN_GENES>>& traits_to_genes,
                                                                                  const std::vector<std::string>& index_to_genes,
                                                                                  const std::unordered_multimap<std::string, std::string>& mapping_genes_to_proteins,
                                                                                  const std::unordered_map<std::string, int>& proteins_to_index,
                                                                                  const std::unordered_multimap<std::string, std::string>& adjacency_list_proteins);

std::unordered_map<std::string, std::bitset<NUM_PHEGEN_PROTEOFORMS>> convertProteinSets(const std::unordered_map<std::string, std::bitset<NUM_PHEGEN_PROTEINS>>& traits_to_proteins,
                                                                                        const std::vector<std::string>& index_to_proteins,
                                                                                        const std::unordered_multimap<std::string, std::string>& mapping_proteins_to_proteoforms,
                                                                                        const std::unordered_map<std::string, int>& proteoforms_to_index,
                                                                                        const std::unordered_multimap<std::string, std::string>& adjacency_list_proteoforms);

template <size_t total_num_entities>
void printMembers(std::ostream& output, const std::bitset<total_num_entities>& entity_set, const std::vector<std::string>& index_to_entities) {
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
                                                                                            const float& min_all_modified_ratio,
                                                                                            const float& min_overlap_modified_ratio) {
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

std::unordered_map<std::string, std::string> createTraitNames(const std::unordered_map<std::string, std::bitset<NUM_PHEGEN_GENES>>& traits_to_genes);

struct reactome_networks {
   std::unordered_multimap<std::string, std::string> adjacency_list_proteins;
   std::unordered_multimap<std::string, std::string> adjacency_list_proteoforms;
};

reactome_networks loadReactomeNetworks(std::string_view path_file_gene_search,
                                       std::string_view path_file_protein_search,
                                       std::string_view path_file_proteoform_search);

std::unordered_multimap<std::string, std::string> loadProteinsAdjacencyList(std::string_view search_file_path);

std::unordered_multimap<std::string, std::string> loadProteoformsAdjacencyList(std::string_view search_file_path);

pathway::entities_bimap deductProteinsFromGenes(std::string_view path_file_mapping_proteins_genes,
                                       const std::unordered_map<std::string, int>& genes_to_index,
                                       const std::unordered_multimap<std::string, std::string>& genes_to_proteins);

pathway::entities_bimap deductProteoformsFromProteins(const std::unordered_multimap<std::string, std::string>& proteins_to_proteoforms, const std::unordered_map<std::string, int>& proteins_to_index);

#endif /* OVERLAP_H_ */