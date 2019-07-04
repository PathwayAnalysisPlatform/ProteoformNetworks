#include "phegeni.hpp"

// Creates bimap for genes and traits from PheGenI data file
load_phegeni_entities_result loadPheGenIEntities(std::string_view path_file_PheGenI,
	const double& max_p_value,
	const umsi& reactome_genes_to_index) {
	const auto [index_to_genes, index_to_traits] = loadPheGenIEntitiesIndexToEntitites(path_file_PheGenI.data(), max_p_value, reactome_genes_to_index);
	const auto genes = createBimap(index_to_genes);
	const auto traits = createBimap(index_to_traits);

	std::cerr << "PHEGEN genes: " << index_to_genes.size() << " = " << genes.entities.size() << "\n";
	std::cerr << "PHEGEN traits: " << index_to_traits.size() << " = " << traits.entities.size() << "\n";

	return { genes,  traits };
}

// Creates vector of genes and traits from PheGenI data file
// entities_column argument starts counting at 0
load_phegeni_entities_index_to_entitites_result loadPheGenIEntitiesIndexToEntitites(const std::string& path_file_PheGenI,
	const double& max_p_value,
	const umsi& reactome_genes_to_index) {
	std::ifstream file_PheGenI(path_file_PheGenI);
	std::string line, field, trait, gene, gene2;
	std::string p_value_str;
	long double p_value;
	uss temp_gene_set, temp_trait_set;
	vs index_to_genes, index_to_traits;

	if (!file_PheGenI.is_open()) {
		throw std::runtime_error("Cannot open " + path_file_PheGenI);
	}

	getline(file_PheGenI, line);                  // Read header line
	while (getline(file_PheGenI, field, '\t')) {  // Read #
		getline(file_PheGenI, trait, '\t');        // Read Trait
		getline(file_PheGenI, field, '\t');        // Read SNP rs
		getline(file_PheGenI, field, '\t');        // Read Context
		getline(file_PheGenI, gene, '\t');         //	Gene
		getline(file_PheGenI, field, '\t');        //	Gene ID
		getline(file_PheGenI, gene2, '\t');        //	Gene 2
		getline(file_PheGenI, field, '\t');        //	Gene ID 2
		getline(file_PheGenI, field, '\t');        // Read Chromosome
		getline(file_PheGenI, field, '\t');        // Read Location
		getline(file_PheGenI, p_value_str, '\t');  // Read P-Value
		getline(file_PheGenI, line);               // Skip header line leftoever: Source,	PubMed,	Analysis ID,	Study ID,	Study Name

		try {
			p_value = stold(p_value_str);
			if (p_value <= GENOME_WIDE_SIGNIFICANCE) {
				temp_trait_set.insert(trait);
				if (reactome_genes_to_index.find(gene) != reactome_genes_to_index.end())
					temp_gene_set.insert(gene);
				if (reactome_genes_to_index.find(gene2) != reactome_genes_to_index.end())
					temp_gene_set.insert(gene2);
			}
		}
		catch (const std::exception& ex) {
			std::cerr << "Error converting: **" << p_value_str << "**\n";
		}
	}

	//Create the final data structures
	index_to_genes = convert_uss_to_vs(temp_gene_set);
	index_to_traits = convert_uss_to_vs(temp_trait_set);

	return { index_to_genes, index_to_traits };
}