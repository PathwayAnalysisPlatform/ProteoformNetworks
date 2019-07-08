#include "overlap.hpp"

using namespace std;

void doOverlapAnalysis(
	std::string_view path_file_PheGenI, 
	std::string_view path_file_genes, 
	std::string_view path_file_mapping_proteins_to_genes, 
	std::string_view path_file_proteoform_search) {

	cout << "Load Reactome genes\n";
	const bimap reactome_genes = createBimap(path_file_genes, true);
	std::cerr << "Number of genes in Reactome: " << reactome_genes.entities.size() << "\n";

	cout << "Loading PheGen genes\n\n";
	const auto [phegeni_genes, phegeni_traits] = loadPheGenIEntities(path_file_PheGenI, reactome_genes);

	cout << "Loading trait gene associations...\n\n";
	const auto mapping_traits_genes = loadPheGenISets(path_file_PheGenI, reactome_genes, phegeni_genes, phegeni_traits);

	// Load mapping from genes to proteins
	cout << "Loading mapping from genes to proteins...\n\n";
	const auto [gene_to_proteins, protein_to_genes] = loadMappingGenesProteins(path_file_mapping_proteins_to_genes);

	// Load mapping from proteins to proteoforms
	cout << "Loading mapping from proteins to proteoforms...\n\n";
	const auto [protein_to_proteoforms, proteoform_to_proteins] = loadMappingProteinsProteoforms(path_file_proteoform_search);

	// Convert gene sets into proteoform sets

	// Calculate overlap score between all trait gene set pairs

	// Calculate overlap score between all trait proteoform set pairs
}

void printMembers(std::ostream& output, const uss& members) {
	int printed = 0;
	output << "[";
	for (auto it = members.begin(); it != members.end(); it++) {
		output << "\"" << *it << "\"";
		if (next(it) != members.end()) {
			output << ",";
		}
	}
	output << "]";
}