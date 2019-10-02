#include <string>
#include <cstdio>
#include <iostream>
#include "overlap_analysis.hpp"

// Input files
const std::string path_file_PheGenI = "../../../resources/PheGenI/PheGenI_Association_genome_wide_significant.txt";
const std::string path_file_reactome_genes = "../../../resources/Reactome/v70/Genes/all_genes_v70.csv";
const std::string path_file_reactome_proteins = "../../../resources/Reactome/v70/Proteins/all_proteins_v70.csv";
const std::string path_file_reactome_proteoforms = "../../../resources/Reactome/v70/Proteoforms/proteoformSearch.tsv";  // Uses this file so that PathwayMatcher has sorted the PTMs in each proteoform.
const std::string path_file_mapping_proteins_to_genes = "../../../resources/UniProt/mapping_proteins_to_genes_v70.tab";
const std::string path_file_mapping_proteins_to_proteoforms = "../../../resources/Reactome/v70/Proteoforms/proteoformSearch.tsv";
const std::string path_file_gene_interactions = "../../../resources/Reactome/v70/Genes/geneInternalEdges.tsv";
const std::string path_file_protein_interactions = "../../../resources/Reactome/v70/Proteins/proteinEdges.tsv";
const std::string path_file_proteoform_interactions = "../../../resources/Reactome/v70/Proteoforms/proteoformEdges.tsv";
const std::string path_scores = "../../../reports/";
const std::string path_modules = "../../../reports/modules/";

int main() try {
//	auto out = freopen("out.txt", "w", stdout);
//	auto err = freopen("err.txt", "w", stderr);

    doOverlapAnalysis(path_file_PheGenI,
                      path_file_reactome_genes,
                      path_file_reactome_proteins,
                      path_file_reactome_proteoforms,
                      path_file_mapping_proteins_to_genes,
                      path_file_mapping_proteins_to_proteoforms,
                      path_file_gene_interactions,
                      path_file_protein_interactions,
                      path_file_proteoform_interactions,
                      path_scores,
                      path_modules);

}
catch (const std::exception &ex) {
    std::cout << ex.what() << "\n";
}
