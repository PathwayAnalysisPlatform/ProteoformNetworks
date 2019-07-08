#include "phegeni.hpp"

using namespace std;

// Creates bimap for genes and traits from PheGenI data file
load_phegeni_genes_and_traits_result loadPheGenIGenesAndTraits(
	string_view path_file_phegeni,
	const bimap& reactome_genes) {

	ifstream file_PheGenI(path_file_phegeni.data());
	string line, field, trait, gene, gene2;
	string p_value_str;
	long double p_value;
	uss temp_gene_set, temp_trait_set;

	std::cerr << "Loaded " << reactome_genes.entities.size() << " = " << reactome_genes.indexes.size() << " reactome genes.\n\n";

	if (!file_PheGenI.is_open()) {
		throw runtime_error("Cannot open path_file_phegeni at " __FUNCTION__);
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

	cerr << "PHEGEN genes: " << phegeni_genes.indexes.size() << " = " << phegeni_genes.entities.size() << "\n";
	cerr << "PHEGEN traits: " << phegeni_traits.indexes.size() << " = " << phegeni_traits.entities.size() << "\n";

	return { phegeni_genes,  phegeni_traits };
}

load_phegeni_sets_result loadPheGenISets(
	string_view path_file_phegeni,
	const bimap& reactome_genes,
	const bimap& phegeni_genes,
	const bimap& phegeni_traits) {
	ifstream file_phegen(path_file_phegeni.data());
	string line, field, trait, gene, gene2;
	string p_value_str;
	long double p_value;
	phegeni_trait_to_genes traits_to_genes;
	phegeni_gene_to_traits genes_to_traits;

	if (!file_phegen.is_open()) {
		throw runtime_error("Cannot open path_file_phegeni at " __FUNCTION__);
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
			traits_to_genes[trait].set(phegeni_genes.indexes.at(gene));
			genes_to_traits[gene].set(phegeni_traits.indexes.at(trait));
		}
		if (reactome_genes.indexes.find(gene2) != reactome_genes.indexes.end()) {
			traits_to_genes[trait].set(phegeni_genes.indexes.at(gene2));
			genes_to_traits[gene2].set(phegeni_traits.indexes.at(trait));
		}
	}

	cerr << "Number of traits with gene members as bitset: " << traits_to_genes.size() << "\n";
	cerr << "Number of genes with traits they belong as bitset: " << genes_to_traits.size() << "\n";

	return { traits_to_genes, genes_to_traits };
}

phegeni_trait_to_proteins convertGeneSets(
	const phegeni_trait_to_genes& traits_to_genes,
	const bimap& phegeni_genes,
	const ummss& mapping_genes_to_proteins,
	const bimap& proteins,
	const ummss& adjacency_list_proteins) {
	cout << "Converting gene to protein sets.\n";
	return convertSets<PHEGENI_GENES, PHEGENI_PROTEINS>(traits_to_genes, phegeni_genes.entities, mapping_genes_to_proteins, proteins.indexes, adjacency_list_proteins);
}

phegeni_trait_to_proteoforms convertProteinSets(const phegeni_trait_to_proteins& traits_to_proteins,
	const bimap& proteins,
	const ummss& mapping_proteins_to_proteoforms,
	const bimap& proteoforms,
	const ummss& adjacency_list_proteoforms) {
	cout << "Converting protein to proteoform sets.\n";
	return convertSets<PHEGENI_PROTEINS, PHEGENI_PROTEOFORMS>(traits_to_proteins, proteins.entities, mapping_proteins_to_proteoforms, proteoforms.indexes, adjacency_list_proteoforms);
}

//uss getGeneStrings(const bitset<PHEGENI_GENES>& gene_set, const bimap& genes) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (gene_set.test(I)) {
//			result.insert(genes.entities[I]);
//		}
//	}
//	return result;
//}
//
//uss getProteinStrings(const bitset<PHEGENI_PROTEINS>& protein_set, const bimap& proteins) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (protein_set.test(I)) {
//			result.insert(proteins.entities[I]);
//		}
//	}
//	return result;
//}
//
//uss getProteoformStrings(const bitset<PHEGENI_PROTEOFORMS>& proteoform_set, const bimap& proteoforms) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (proteoform_set.test(I)) {
//			result.insert(proteoforms.entities[I]);
//		}
//	}
//	return result;
//}
//

umss createTraitNames(const um<string, bitset<PHEGENI_GENES>>& traits_to_genes) {
	umss sets_to_names;
	for (const auto& trait_entry : traits_to_genes) {
		sets_to_names.emplace(trait_entry.first, trait_entry.first);
	}
	return sets_to_names;
}