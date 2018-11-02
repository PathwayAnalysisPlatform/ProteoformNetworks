#ifndef ARTEFACTUAL_OVERLAP_H_
#define ARTEFACTUAL_OVERLAP_H_

#include <set>
#include <map>
#include <bitset>
#include <string>
#include <vector>

void doArtefactualOverlapAnalysis(const std::string& report_path);

std::set<std::pair<std::string, std::string>> findDiseaseModulePairsWithArtefactualOverlap(const std::string& path_file_gene_art_pairs,
                                                                                           const std::string& path_file_protein_art_pairs);

void doArtefactualOverlapAnalysis(const std::string& path_file_gene_search,
                                  const std::string& path_file_protein_search,
                                  const std::string& path_file_proteoform_search,
                                  const std::string& report_path,
                                  const std::string& path_file_gene_art_pairs,
                                  const std::string& path_file_protein_art_pairs);

#endif /* ARTEFACTUAL_OVERLAP_H_ */
