#include <string>
#include <cstdio>
#include <iostream>

#include "overlap.hpp"

// Input files
const std::string path_file_reactome_genes = "resources/Reactome/v69/genes.csv";
const std::string path_file_PheGenI = "resources/PheGenI/PheGenI_Association_full.tab";
//const std::string path_file_protein_search = "resources/Reactome/v69/protein_search.tsv";
//const std::string path_file_proteoform_search = "resources/Reactome/v69/proteoform_search.tsv";
//const std::string path_file_mapping_proteins_to_genes = "resources/UniProt/proteins_to_genes.tab";

int main() try {

	auto out = freopen("out.txt", "w", stdout);
	auto err = freopen("err.txt", "w", stderr);

	doOverlapAnalysis(path_file_PheGenI, path_file_reactome_genes);
}
catch (const std::exception& ex) {
	std::cout << ex.what() << "\n";
}