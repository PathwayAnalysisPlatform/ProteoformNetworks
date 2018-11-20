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

const std::regex RGX_ACCESSION_DELIMITER{"[;-]"};

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
const int MAX_OVERLAP_SIZE = 10;

const int MIN_SET_SIZE = 1;
const int MAX_SET_SIZE = 20;

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

struct Entities_bimap {
   std::vector<std::string> index_to_entities;
   std::unordered_map<std::string, int> entities_to_index;
};

struct Frequencies {
   std::map<std::string, int> modifications;
   std::map<std::string, int> proteins;
   std::map<std::string, int> proteoforms;
};

struct index_to_entitites_phegen_result {
   std::vector<std::string> index_to_genes;
   std::vector<std::string> index_to_traits;
};

struct load_entitites_phegen_result {
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

Entities_bimap loadEntities(const std::string& entities_file_path);

std::unordered_map<std::string, std::string> loadPathwayNames(const std::string& path_protein_search_file);

std::map<std::string, int> fillMap(const std::vector<std::string>& index_to_entities);

std::unordered_map<std::string, std::bitset<NUM_GENES>> loadGeneSets(const std::string& file_path,
                                                                     const std::unordered_map<std::string, int>& entities_to_index,
                                                                     bool pathways);
std::unordered_map<std::string, std::bitset<NUM_PROTEINS>> loadProteinSets(const std::string& file_path,
                                                                           const std::unordered_map<std::string, int>& entities_to_index,
                                                                           bool pathways);
std::unordered_map<std::string, std::bitset<NUM_PROTEOFORMS>> loadProteoformSets(const std::string& file_path,
                                                                                 const std::unordered_map<std::string, int>& entities_to_index,
                                                                                 bool pathways);

std::map<std::pair<std::string, std::string>, std::bitset<NUM_GENES>> findOverlappingGeneSets(const std::unordered_map<std::string, std::bitset<NUM_GENES>>& sets_to_members,
                                                                                              const int& min_overlap, const int& max_overlap,
                                                                                              const int& min_set_size, const int& max_set_size);
std::map<std::pair<std::string, std::string>, std::bitset<NUM_PROTEINS>> findOverlappingProteinSets(const std::unordered_map<std::string, std::bitset<NUM_PROTEINS>>& sets_to_members,
                                                                                                    const int& min_overlap, const int& max_overlap,
                                                                                                    const int& min_set_size, const int& max_set_size);

std::map<std::pair<std::string, std::string>, std::bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const std::unordered_map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members);
std::map<std::pair<std::string, std::string>, std::bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const std::unordered_map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                                                          const int& min_overlap, const int& max_overlap,
                                                                                                          const int& min_set_size, const int& max_set_size);
std::map<std::pair<std::string, std::string>, std::bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const std::unordered_map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                                                          const int& min_overlap, const int& max_overlap,
                                                                                                          const int& min_set_size, const int& max_set_size,
                                                                                                          const std::bitset<NUM_PROTEOFORMS>& modified_proteoforms,
                                                                                                          const float& min_all_modified_ratio,
                                                                                                          const float& min_overlap_modified_ratio);

void printGeneMembers(std::ofstream& output,
                      const std::bitset<NUM_GENES>& gene_set,
                      const std::vector<std::string>& index_to_entities);
void printProteinMembers(std::ofstream& output,
                         const std::bitset<NUM_PROTEINS>& protein_set,
                         const std::vector<std::string>& index_to_entities);
void printProteoformMembers(std::ofstream& output,
                            const std::bitset<NUM_PROTEOFORMS>& proteoform_set,
                            const std::vector<std::string>& index_to_entities);

std::set<std::string> getGeneStrings(const std::bitset<NUM_GENES>& gene_set,
                                     const std::vector<std::string>& index_to_genes);
std::set<std::string> getProteinStrings(const std::bitset<NUM_PROTEINS>& protein_set,
                                        const std::vector<std::string>& index_to_proteins);
std::set<std::string> getProteoformStrings(const std::bitset<NUM_PROTEOFORMS>& proteoform_set,
                                           const std::vector<std::string>& index_to_proteoforms);

std::string getAccession(std::string proteoform);

load_entitites_phegen_result loadEntitiesPheGen(const std::string& path_file_PheGenI_full,
                                                const double& max_p_value);

load_trait_gene_sets_result loadTraitGeneSets(const std::string& path_file_phegen,
                                              const double& max_p_value,
                                              const std::vector<std::string>& index_to_genes,
                                              const std::vector<std::string>& index_to_traits,
                                              const std::unordered_map<std::string, int>& genes_to_index,
                                              const std::unordered_map<std::string, int>& traits_to_index);

load_mapping_result loadMapping(const std::string& path_file_mapping);

std::unordered_map<std::string, std::bitset<NUM_PHEGEN_PROTEINS>> convertGeneSetsToProteinSets(const std::unordered_map<std::string, std::bitset<NUM_PHEGEN_GENES>>& traits_to_genes,
                                                                                               const std::vector<std::string>& index_to_genes,
                                                                                               const std::unordered_multimap<std::string, std::string>& mapping_genes_to_proteins,
                                                                                               const std::unordered_map<std::string, int>& proteins_to_index,
                                                                                               const std::unordered_multimap<std::string, std::string>& adjacency_list_proteins);

std::vector<std::string> convert(const std::unordered_set<std::string>& a_set);

std::unordered_map<std::string, int> getEntitiesToIndex(const std::vector<std::string>& index_to_entities);

#endif /* OVERLAP_H_ */