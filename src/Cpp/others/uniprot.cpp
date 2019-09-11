#include "uniprot.hpp"

// Loads mapping from one column of strings to a second column of strings. If there are multiple values on the second column for one value on the left, then
//they are separated by spaces. Columns are separated by tabs.
load_mapping_genes_proteins_result loadMappingGenesProteins(std::string_view path_file_mapping) {
	ummss gene_to_proteins;
	ummss protein_to_genes;
	std::ifstream file_mapping(path_file_mapping.data());
	std::string gene, protein, leftover;

	if (!file_mapping.is_open()) {
		std::string message = "Error reading file at: ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
	}

	getline(file_mapping, leftover);  // Discard the header line.
	while (getline(file_mapping, protein, '\t')) {
		char c;

		while (file_mapping.get(c)) {
			if (c == ' ' || c == '\t' || c == '\n') {
				gene_to_proteins.emplace(gene, protein);
				protein_to_genes.emplace(protein, gene);
				gene = "";
				if (c == '\t') {
					getline(file_mapping, leftover);
				}
				if (c != ' ') {
					break;
				}
			}
			else {
				gene = gene.append(1, c);
			}
		}
	}

	std::cerr << "Mapping loaded: " << gene_to_proteins.size() << " = " << protein_to_genes.size() << "\n";
	return { gene_to_proteins, protein_to_genes };
}