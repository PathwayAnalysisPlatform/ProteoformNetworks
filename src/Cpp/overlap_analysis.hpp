#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <iostream>
#include <ctime>
#include <sys/stat.h>
#include <utility>

#include "phegeni.hpp"
#include "scores.hpp"
#include "bimap_str_int.hpp"
#include "others/uniprot.hpp"
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

struct get_entities_result {
    bimap_str_int genes;
    bimap_str_int proteins;
    bimap_str_int proteoforms;
};

get_entities_result get_entities(std::string_view path_file_reactome_genes,
                                 std::string_view path_file_reactome_proteins,
                                 std::string_view path_file_reactome_proteoforms);

inline bool file_exists(const std::string &name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

struct get_modules_result {
    modules gene_modules;
    modules protein_modules;
    modules proteoform_modules;
};

// Create or read module files at the three levels: all in one, and single module files.
get_modules_result get_or_create_modules(std::string path_modules,
                                         std::string_view path_file_phegeni,
                                         std::string_view path_file_gene_interactions,
                                         std::string_view path_file_mapping_proteins_to_genes,
                                         std::string_view path_file_protein_interactions,
                                         std::string_view path_file_mapping_proteins_to_proteoforms,
                                         std::string_view path_file_proteoform_interactions,
                                         const bimap_str_int &genes,
                                         const bimap_str_int &proteins,
                                         const bimap_str_int &proteoforms,
                                         const bimap_str_int &traits);

struct score_maps {
    pair_map<double> gene_scores;
    pair_map<double> protein_scores;
    pair_map<double> proteoform_scores;
};

score_maps get_scores(std::string path_scores,
                      std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>)> scoring,
                      std::string label,
                      const modules &gene_modules,
                      const modules &protein_modules,
                      const modules &proteoform_modules,
                      const bimap_str_int &traits);

// Calculates differences in module overlap between gene and proteoform level networks_lib
// It uses multiple scoring functions to calculate the overlap score between each pair of modules.
// Creates a report for each scoring function. The report is a csv file with one row each module pair.
void doOverlapAnalysis(
        std::string_view path_file_phegeni,
        std::string_view path_file_reactome_genes,
        std::string_view path_file_reactome_proteins,
        std::string_view path_file_reactome_proteoforms,
        std::string_view path_file_mapping_proteins_to_genes,
        std::string_view path_file_mapping_proteins_to_proteoforms,
        std::string_view path_file_gene_interactions,
        std::string_view path_file_protein_interactions,
        std::string_view path_file_proteoform_interactions,
        std::string path_reports,
        std::string_view path_modules);

void report_module_size_variation(std::string_view path_reports, const modules &gene_modules,
                                  const modules &protein_modules, const modules &proteoform_modules,
                                  const bimap_str_int &traits);

// Check if there are pairs of modules sharing nodes at one level but not in another
void report_node_overlap_reduction_examples(std::string path_scores, std::string label, const score_maps &scores,
                                            const bimap_str_int &traits);

// Check if there are pairs of modules which are have variations in level, for the number of edges connecting them
void report_connecting_edges_variation_examples(std::string path_reports, const score_maps scores);

void report_overlap_only_ptms(std::string string, const score_maps maps, const bimap_str_int anInt);

// Version with set size and overlap size limits
template<size_t total_num_entities>
std::map<std::pair<std::string, std::string>, std::bitset<total_num_entities>>
findOverlappingPairs(const um<std::string, std::bitset<total_num_entities>> &sets_to_members,
                     const int &min_overlap, const int &max_overlap,
                     const int &min_set_size, const int &max_set_size) {
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
    std::cerr << "tardamos " << (t1 - static_cast<double>(t0)) / CLOCKS_PER_SEC << "\n";
    return result;
}

// This version includes modified ratio of the proteoform members of the set, and another ratio for the overlapping proteoforms.
template<size_t S>
std::map<std::pair<std::string, std::string>, std::bitset<S>>
findOverlappingProteoformSets(const std::unordered_map<std::string, std::bitset<S>> &sets_to_members,
                              const int &min_overlap, const int &max_overlap,
                              const int &min_set_size, const int &max_set_size,
                              const std::bitset<S> &modified_proteoforms,
                              const double &min_all_modified_ratio,
                              const double &min_overlap_modified_ratio) {
    std::vector<typename std::unordered_map<std::string, std::bitset<S>>::const_iterator> nav;
    for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
        int set_size = it->second.count();
        if (min_set_size <= set_size && set_size <= max_set_size) {
            //  cerr << "MIN_SET_SIZE: " << min_set_size << " SET_SIZE: " << set_size << " MAX_SET_SIZE: " << max_set_size << "\n";
            float percentage =
                    static_cast<float>((modified_proteoforms & it->second).count()) / static_cast<float>(set_size);
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
                float percentage = static_cast<float>((modified_proteoforms & overlap).count()) /
                                   static_cast<float>(overlap.count());
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

template<size_t total_num_entities>
void
printMembers(std::ostream &output, const std::bitset<total_num_entities> &entity_set, const bimap_str_int &entities) {
    int printed = 0;
    int total = entity_set.count();
    output << "[";
    for (int I = 0; I < total_num_entities; I++) {
        if (entity_set.test(I)) {
            output << "\"" << entities.int_to_str[I] << "\"";
            printed++;
            if (printed != total) {
                output << ",";
            }
        }
    }
    output << "]";
}

void printMembers(std::ostream &output, const uss &members);

template<size_t total_num_entities>
uss getStringSetFromEntityBitset(const std::bitset<total_num_entities> &entity_set, const bimap_str_int &entities) {
    uss result;
    for (int I = 0; I < total_num_entities; I++) {
        if (entity_set.test(I)) {
            result.insert(entities.int_to_str[I]);
        }
    }
    return result;
}

void report_overlaps_with_ptms();

#endif /* OVERLAP_H_ */