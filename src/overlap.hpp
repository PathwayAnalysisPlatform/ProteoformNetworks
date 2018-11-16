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
const size_t NUM_PHEGEN_PROTEINS = 20408;

const int MIN_OVERLAP_SIZE = 1;
const int MAX_OVERLAP_SIZE = 10;

const int MIN_PATHWAY_SIZE = 1;
const int MAX_PATHWAY_SIZE = 20;

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

struct Entities_bimap {
   std::vector<std::string> index_to_entities;
   std::map<std::string, int> entities_to_index;
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
   std::map<std::string, int> genes_to_index;
   std::map<std::string, int> traits_to_index;
};

struct load_trait_gene_sets_result {
   std::map<std::string, std::bitset<NUM_PHEGEN_TRAITS>> genes_to_sets;
   std::map<std::string, std::bitset<NUM_PHEGEN_GENES>> sets_to_genes;
};

Entities_bimap loadEntities(const std::string& entities_file_path);

std::map<std::string, std::string> loadPathwayNames(const std::string& path_protein_search_file);

std::map<std::string, int> fillMap(const std::vector<std::string>& index_to_entities);

std::map<std::string, std::bitset<NUM_GENES>> loadPathwaysGeneMembers(const std::string& file_path, const std::map<std::string, int>& entities_to_index);
std::map<std::string, std::bitset<NUM_PROTEINS>> loadPathwaysProteinMembers(const std::string& file_path, const std::map<std::string, int>& entities_to_index);
std::map<std::string, std::bitset<NUM_PROTEOFORMS>> loadPathwaysProteoformMembers(const std::string& file_path, const std::map<std::string, int>& entities_to_index);

std::map<std::pair<std::string, std::string>, std::bitset<NUM_GENES>> findOverlappingGeneSets(const std::map<std::string, std::bitset<NUM_GENES>>& sets_to_members,
                                                                                              const int& min_overlap, const int& max_overlap,
                                                                                              const int& min_set_size, const int& max_set_size);
std::map<std::pair<std::string, std::string>, std::bitset<NUM_PROTEINS>> findOverlappingProteinSets(const std::map<std::string, std::bitset<NUM_PROTEINS>>& sets_to_members,
                                                                                                    const int& min_overlap, const int& max_overlap,
                                                                                                    const int& min_set_size, const int& max_set_size);

std::map<std::pair<std::string, std::string>, std::bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const std::map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members);
std::map<std::pair<std::string, std::string>, std::bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const std::map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                                                          const int& min_overlap, const int& max_overlap,
                                                                                                          const int& min_set_size, const int& max_set_size);
std::map<std::pair<std::string, std::string>, std::bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const std::map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                                                          const int& min_overlap, const int& max_overlap,
                                                                                                          const int& min_set_size, const int& max_set_size,
                                                                                                          const std::bitset<NUM_PROTEOFORMS>& modified_proteoforms,
                                                                                                          const float& min_all_modified_ratio,
                                                                                                          const float& min_overlap_modified_ratio);

void printGeneMembers(std::ofstream& output, const std::bitset<NUM_GENES>& gene_set, const std::vector<std::string>& index_to_entities);
void printProteinMembers(std::ofstream& output, const std::bitset<NUM_PROTEINS>& protein_set, const std::vector<std::string>& index_to_entities);
void printProteoformMembers(std::ofstream& output, const std::bitset<NUM_PROTEOFORMS>& proteoform_set, const std::vector<std::string>& index_to_entities);

std::set<std::string> getGeneStrings(const std::bitset<NUM_GENES>& gene_set, const std::vector<std::string>& index_to_genes);
std::set<std::string> getProteinStrings(const std::bitset<NUM_PROTEINS>& protein_set, const std::vector<std::string>& index_to_proteins);
std::set<std::string> getProteoformStrings(const std::bitset<NUM_PROTEOFORMS>& proteoform_set, const std::vector<std::string>& index_to_proteoforms);

std::string getAccession(std::string proteoform);

load_entitites_phegen_result loadEntitiesPheGen(const std::string& path_file_PheGenI_full, const double& max_p_value);

load_trait_gene_sets_result loadTraitGeneSets(const std::string& path_file_phegen,
                                              const double& max_p_value,
                                              const std::vector<std::string>& index_to_genes,
                                              const std::vector<std::string>& index_to_traits,
                                              const std::map<std::string, int>& genes_to_index,
                                              const std::map<std::string, int>& traits_to_index);

std::multimap<std::string, std::string> loadMapping(const std::string& path_file_mapping);

template <size_t total_num_proteins>
std::map<std::string, std::bitset<total_num_proteins>> convertGeneSetsToProteinSets(const std::map<std::string, std::bitset<NUM_PHEGEN_GENES>>& traits_to_genes,
                                                                                    const std::multimap<std::string, std::string>& proteins_to_genes);

std::vector<std::string> convert(const std::unordered_set<std::string>& a_set);

#endif /* OVERLAP_H_ */