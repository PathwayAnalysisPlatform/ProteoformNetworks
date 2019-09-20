#include <string>
#include <cstdio>
#include <iostream>

#include "overlap_analysis.hpp"

// Input files
const std::string path_file_PheGenI = "../../../resources/PheGenI/PheGenI_Association_genome_wide_significant.txt";
const std::string path_file_reactome_genes = "../../../resources/Reactome/v70/Genes/all_genes_v70.csv";
const std::string path_file_reactome_proteins = "../../../resources/Reactome/v70/Proteins/all_proteins_v70.csv";
const std::string path_file_reactome_proteoforms = "../../../resources/Reactome/v70/Proteoforms/all_proteoforms_v70.csv";
const std::string path_file_mapping_proteins_to_genes = "../../../resources/UniProt/mapping_protein_to_genes.tab";
const std::string path_file_mapping_proteins_to_proteoforms = "../../../resources/Reactome/v70/mapping_protein_to_genes.tab";
const std::string path_file_protein_interactions = "../../../resources/Reactome/v70/Proteins/proteinInternalEdges.tsv";
const std::string path_file_proteoform_interactions = "../../../resources/Reactome/v70/Proteoforms/proteoformInternalEdges.tsv";
const std::string path_scores = "../../../reports/";

int main() try {
//	auto out = freopen("out.txt", "w", stdout);
//	auto err = freopen("err.txt", "w", stderr);

    doOverlapAnalysis(path_file_PheGenI,
                      path_file_reactome_genes,
                      path_file_reactome_proteins,
                      path_file_reactome_proteoforms,
                      path_file_mapping_proteins_to_genes,
                      path_file_mapping_proteins_to_proteoforms,
                      path_file_protein_interactions,
                      path_file_proteoform_interactions,
                      path_scores);
}
catch (const std::exception &ex) {
    std::cout << ex.what() << "\n";
}
