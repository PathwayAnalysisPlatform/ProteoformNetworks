#include "overlap.hpp"

using namespace std;

void doOverlapAnalysis(std::string_view path_file_PheGenI, std::string_view path_file_genes) {

	cout << "Load Reactome genes\n";
	const bimap reactome_genes = createBimap(path_file_genes, true);

	cout << "Loading PheGen genes\n";
	const auto [phegeni_genes, phegeni_traits] = loadPheGenIEntities(path_file_PheGenI, reactome_genes, GENOME_WIDE_SIGNIFICANCE);

	std::cout << "Number of genes in Reactome: " << reactome_genes.entities.size() << "\n";
	std::cout << "Number of genes in PheGenI: " << phegeni_genes.entities.size() << "\n";
	std::cout << "Number of traits in PheGenI: " << phegeni_traits.entities.size() << "\n";

	const auto mapping_traits_genes = loadPheGenISets(path_file_PheGenI, reactome_genes, GENOME_WIDE_SIGNIFICANCE, phegeni_genes, phegeni_traits);

	std::cout << "Number of traits with gene members as bitset: " << mapping_traits_genes.trait_to_genes.size() << "\n";
	std::cout << "Number of genes with traits they belong as bitset: " << mapping_traits_genes.gene_to_traits.size() << "\n";
}