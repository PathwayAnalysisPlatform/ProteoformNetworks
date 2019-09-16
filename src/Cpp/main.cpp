#include <string>
#include <cstdio>
#include <iostream>

#include "overlap_analysis.hpp"

#include <windows.h>

std::string GetExeFileName() {
    char buffer[MAX_PATH];
    GetModuleFileName(NULL, buffer, MAX_PATH);
    return std::string(buffer);
}

std::string GetExePath() {
    std::string f = GetExeFileName();
    return f.substr(0, f.find_last_of("\\/"));
}

// Input files
const std::string path_file_PheGenI = "../../../resources/PheGenI/PheGenI_Association_genome_wide_significant.txt";
const std::string path_file_reactome_genes = "../../../resources/Reactome/v70/Genes/all_genes_v70.csv";
const std::string path_file_reactome_proteins = "../../../resources/Reactome/v70/Proteins/all_proteins_v70.csv";
const std::string path_file_mapping_proteins_to_genes = "../../../resources/UniProt/mapping_protein_to_genes.tab";
const std::string path_file_protein_edges = "../../../resources/Reactome/v70/Proteins/proteinInternalEdges.tsv";
const std::string path_file_proteoform_edges = "../../../resources/Reactome/v70/Proteoforms/proteoformInternalEdges.tsv";
const std::string path_scores = "../../../reports/";

int main() try {
//	auto out = freopen("out.txt", "w", stdout);
//	auto err = freopen("err.txt", "w", stderr);
    std::cout << GetExePath() << std::endl;

    doOverlapAnalysis(path_file_PheGenI,
                      path_file_reactome_genes,
                      path_file_reactome_proteins,
                      path_file_mapping_proteins_to_genes,
                      path_file_protein_edges,
                      path_file_proteoform_edges,
                      path_scores);
}
catch (const std::exception &ex) {
    std::cout << ex.what() << "\n";
}
