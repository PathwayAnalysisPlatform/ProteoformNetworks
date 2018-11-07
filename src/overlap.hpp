#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <algorithm>
#include <bitset>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

const size_t NUM_GENES = 23970;
const size_t NUM_PROTEINS = 10778;
const size_t NUM_PROTEOFORMS = 13911;

const int MIN_OVERLAP_SIZE = 1;
const int MAX_OVERLAP_SIZE = 10;

const int MIN_PATHWAY_SIZE = 1;
const int MAX_PATHWAY_SIZE = 20;

struct Entities_bimap {
    std::vector<std::string> index_to_entities;
    std::map<std::string, int> entities_to_index;
};

Entities_bimap loadEntities(const std::string& entities_file_path);

std::map<std::string, std::string> loadPathwayNames(const std::string& path_protein_search_file);

std::map<std::string, int> fillMap(const std::vector<std::string>& index_to_entities);

std::map<std::string, std::bitset<NUM_GENES>> loadPathwaysGeneMembers(const std::string& file_path, const std::map<std::string, int>& entities_to_index);
std::map<std::string, std::bitset<NUM_PROTEINS>> loadPathwaysProteinMembers(const std::string& file_path, const std::map<std::string, int>& entities_to_index);
std::map<std::string, std::bitset<NUM_PROTEOFORMS>> loadPathwaysProteoformMembers(const std::string& file_path, const std::map<std::string, int>& entities_to_index);

std::set<std::pair<std::string, std::string>> findOverlappingGeneSets(const std::map<std::string, std::bitset<NUM_GENES>>& sets_to_members,
                                                                      const int& min_overlap, const int& max_overlap,
                                                                      const int& min_set_size, const int& max_set_size);
std::set<std::pair<std::string, std::string>> findOverlappingProteinSets(const std::map<std::string, std::bitset<NUM_PROTEINS>>& sets_to_members,
                                                                         const int& min_overlap, const int& max_overlap,
                                                                         const int& min_set_size, const int& max_set_size);
std::set<std::pair<std::string, std::string>> findOverlappingProteoformSets(const std::map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                            const int& min_overlap, const int& max_overlap,
                                                                            const int& min_set_size, const int& max_set_size);

std::set<std::pair<std::string, std::string>> findOverlappingProteoformSets(const std::map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members);

std::set<std::pair<std::string, std::string>> findOverlappingProteoformSets(const std::map<std::string, std::bitset<NUM_PROTEOFORMS>>& sets_to_members,
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

#endif /* OVERLAP_H_ */