#include <string>
#include <cstdio>
#include <iostream>

#include "overlap.hpp"

// Input files
const std::string path_file_PheGenI = "resources/PheGenI/PheGenI_Association_genome_wide_significant.txt";
const std::string path_file_reactome_genes = "resources/Reactome/v69/genes.csv";
const std::string path_file_mapping_proteins_to_genes = "resources/UniProt/mapping_protein_to_genes.tab";
const std::string path_file_protein_search = "resources/Reactome/v69/protein_search.tsv";
const std::string path_file_proteoform_search = "resources/Reactome/v69/proteoforms_search.tsv";

int main() try {
//	auto out = freopen("out.txt", "w", stdout);
//	auto err = freopen("err.txt", "w", stderr);

	doOverlapAnalysis(path_file_PheGenI, path_file_reactome_genes, path_file_mapping_proteins_to_genes,
	        path_file_protein_search, path_file_proteoform_search);
}
catch (const std::exception& ex) {
	std::cout << ex.what() << "\n";
}
