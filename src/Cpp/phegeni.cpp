#include "phegeni.hpp"

// Creates bimap for genes and traits from PheGenI data file
load_phegeni_entities_result loadPheGenIEntities(
	std::string_view path_file_phegeni,
	const bimap& reactome_genes,
	const double& max_p_value) {

	std::ifstream file_PheGenI(path_file_phegeni.data());
	std::string line, field, trait, gene, gene2;
	std::string p_value_str;
	long double p_value;
	uss temp_gene_set, temp_trait_set;

	if (!file_PheGenI.is_open()) {
		throw std::runtime_error("Cannot open path_file_phegeni at " __FUNCTION__);
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
		getline(file_PheGenI, line);        // Read rest of line

		if (reactome_genes.indexes.find(gene) != reactome_genes.indexes.end()) {
			temp_trait_set.insert(trait);
			temp_gene_set.insert(gene);
		}
		if (reactome_genes.indexes.find(gene2) != reactome_genes.indexes.end()) {
			temp_trait_set.insert(trait);
			temp_gene_set.insert(gene2);
		}
	}
	const auto phegeni_genes = createBimap(convert_uss_to_vs(temp_gene_set));
	const auto phegeni_traits = createBimap(convert_uss_to_vs(temp_trait_set));

	std::cerr << "PHEGEN genes: " << phegeni_genes.indexes.size() << " = " << phegeni_genes.entities.size() << "\n";
	std::cerr << "PHEGEN traits: " << phegeni_traits.indexes.size() << " = " << phegeni_traits.entities.size() << "\n";

	return { phegeni_genes,  phegeni_traits };
}

load_phegeni_sets_result loadPheGenISets(
	std::string_view path_file_phegeni,
	const bimap& reactome_genes,
	const double& max_p_value,
	const bimap& phegeni_genes,
	const bimap& phegeni_traits) {
	std::ifstream file_phegen(path_file_phegeni.data());
	std::string line, field, trait, gene, gene2;
	std::string p_value_str;
	long double p_value;
	phegeni_trait_to_genes trait_to_genes;
	phegeni_gene_to_traits gene_to_traits;

	if (!file_phegen.is_open()) {
		throw std::runtime_error("Cannot open path_file_phegeni at " __FUNCTION__);
	}

	getline(file_phegen, line);                  // Read header line
	while (getline(file_phegen, field, '\t')) {  // Read #
		getline(file_phegen, trait, '\t');        // Read Trait
		getline(file_phegen, field, '\t');        // Read SNP rs
		getline(file_phegen, field, '\t');        // Read Context
		getline(file_phegen, gene, '\t');         //	Gene
		getline(file_phegen, field, '\t');        //	Gene ID
		getline(file_phegen, gene2, '\t');        //	Gene 2
		getline(file_phegen, field, '\t');        //	Gene ID 2
		getline(file_phegen, field, '\t');        // Read Chromosome
		getline(file_phegen, field, '\t');        // Read Location
		getline(file_phegen, p_value_str, '\t');  // Read P-Value
		getline(file_phegen, line);               // Skip header line leftoever: Source,	PubMed,	Analysis ID,	Study ID,	Study Name

		if (reactome_genes.indexes.find(gene) != reactome_genes.indexes.end()) {
			trait_to_genes[trait].set(phegeni_genes.indexes.at(gene));
			gene_to_traits[gene].set(phegeni_traits.indexes.at(trait));
		}
		if (reactome_genes.indexes.find(gene2) != reactome_genes.indexes.end()) {
			trait_to_genes[trait].set(phegeni_genes.indexes.at(gene2));
			gene_to_traits[gene2].set(phegeni_traits.indexes.at(trait));
		}
	}

	return { trait_to_genes, gene_to_traits };
}

//uss getGeneStrings(const std::bitset<PHEGENI_GENES>& gene_set, const bimap& genes) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (gene_set.test(I)) {
//			result.insert(genes.entities[I]);
//		}
//	}
//	return result;
//}
//
//uss getProteinStrings(const std::bitset<PHEGENI_PROTEINS>& protein_set, const bimap& proteins) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (protein_set.test(I)) {
//			result.insert(proteins.entities[I]);
//		}
//	}
//	return result;
//}
//
//uss getProteoformStrings(const std::bitset<PHEGENI_PROTEOFORMS>& proteoform_set, const bimap& proteoforms) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (proteoform_set.test(I)) {
//			result.insert(proteoforms.entities[I]);
//		}
//	}
//	return result;
//}
//
//umss createTraitNames(const um<std::string, std::bitset<PHEGENI_GENES>>& traits_to_genes) {
//	umss sets_to_names;
//	for (const auto& trait_entry : traits_to_genes) {
//		sets_to_names.emplace(trait_entry.first, trait_entry.first);
//	}
//	return sets_to_names;
//}