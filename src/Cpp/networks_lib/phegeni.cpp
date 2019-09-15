#include <cstring>
#include <iostream>
#include "phegeni.hpp"

using namespace std;

// Read the genes and traits of PheGenI data file.
// Creates two bimaps: genes and traits
// The gene list are those genes in the dataset that are also contained in the acceptable gene list.
load_phegeni_genes_and_traits_result loadPheGenIGenesAndTraits(
        string_view path_file_phegeni,
        const bimap_str_int &acceptable_genes) {

    ifstream file_PheGenI(path_file_phegeni.data());
    string line, field, trait, gene, gene2;
    string p_value_str;
    long double p_value;
    uss temp_gene_set, temp_trait_set;

    if (!file_PheGenI.is_open()) {
        std::string message = "Cannot open path_file_phegeni at ";
        std::string function = __FUNCTION__;
        throw runtime_error(message + function);
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

        if (acceptable_genes.str_to_int.find(gene) != acceptable_genes.str_to_int.end()) {
            temp_trait_set.insert(trait);
            temp_gene_set.insert(gene);
        }
        if (acceptable_genes.str_to_int.find(gene2) != acceptable_genes.str_to_int.end()) {
            temp_trait_set.insert(trait);
            temp_gene_set.insert(gene2);
        }
    }
    vs gene_vector = convert_uss_to_vs(temp_gene_set);
    sort(gene_vector.begin(), gene_vector.end());
    const auto phegeni_genes = createBimap(gene_vector);

    vs traits_vector = convert_uss_to_vs(temp_trait_set);
    sort(traits_vector.begin(), traits_vector.end());
    const auto phegeni_traits = createBimap(traits_vector);

    cerr << "PHEGEN genes: " << phegeni_genes.str_to_int.size() << " = " << phegeni_genes.int_to_str.size() << "\n";
    cerr << "PHEGEN traits: " << phegeni_traits.str_to_int.size() << " = " << phegeni_traits.int_to_str.size() << "\n";

    return {phegeni_genes, phegeni_traits};
}

// Creates two mappings to represent the Trait modules with genes as members:
// - trait strings to gene bitsets.
// - gene strings to trait bitsets
// The two mappings define the trait modules (also called, phenotype modules or disease modules)
trait_modules loadPheGenIModules(
        string_view path_file_phegeni,
        const bimap_str_int &acceptable_genes) {
    ifstream file_phegen(path_file_phegeni.data());
    string line, field, trait, gene, gene2;
    string p_value_str;
    long double p_value;

    if (!file_phegen.is_open()) {
        std::string message = "Cannot open path_file_phegeni at ";
        std::string function = __FUNCTION__;
        throw runtime_error(message + function);
    }

    const auto[phegeni_genes, phegeni_traits] = loadPheGenIGenesAndTraits(path_file_phegeni, acceptable_genes);

    // Initialize modules
    trait_modules modules;
    for (const auto &gene : phegeni_genes.int_to_str) {
        modules.genes_to_traits[gene] = base::dynamic_bitset<>(phegeni_traits.int_to_str.size());
    }
    for (const auto &trait : phegeni_traits.int_to_str) {
        modules.traits_to_genes[trait] = base::dynamic_bitset<>(phegeni_genes.int_to_str.size());
    }

    std::cout << "Initialized modules:\n"
              << "Genes to traits: " << modules.genes_to_traits.size()
              << "\t Traits bitsets size: " << modules.genes_to_traits.begin()->second.size()
              << "\nTraits to genes: " << modules.traits_to_genes.size()
              << "\t Genes bitsets size: " << modules.traits_to_genes.begin()->second.size() << "\n";

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

        if (acceptable_genes.str_to_int.find(gene) != acceptable_genes.str_to_int.end()) {
            modules.traits_to_genes[trait][phegeni_genes.str_to_int.at(gene)].set();
            modules.genes_to_traits[gene][phegeni_traits.str_to_int.at(trait)].set();
        }
        if (acceptable_genes.str_to_int.find(gene2) != acceptable_genes.str_to_int.end()) {
            modules.traits_to_genes[trait][phegeni_genes.str_to_int.at(gene2)].set();
            modules.genes_to_traits[gene2][phegeni_traits.str_to_int.at(trait)].set();
        }
    }

    cerr << "Number of traits with gene members as bitset: " << modules.traits_to_genes.size() << "\n";
    cerr << "Number of genes with traits they belong as bitset: " << modules.genes_to_traits.size() << "\n";

    return modules;
}

umsb convertGeneSets(
        const umsb &traits_to_genes,
        const bimap_str_int &phegeni_genes,
        const ummss &mapping_genes_to_proteins,
        const bimap_str_int &proteins,
        const ummss &adjacency_list_proteins) {
    cout << "Converting gene to protein sets.\n";
    return convertSets(traits_to_genes, phegeni_genes.int_to_str,
                       mapping_genes_to_proteins, proteins.str_to_int,
                       adjacency_list_proteins);
}

umsb convertProteinSets(const umsb &traits_to_proteins,
                        const bimap_str_int &proteins,
                        const ummss &mapping_proteins_to_proteoforms,
                        const bimap_str_int &proteoforms,
                        const ummss &adjacency_list_proteoforms) {
    cout << "Converting protein to proteoform sets.\n";
    return convertSets(traits_to_proteins, proteins.int_to_str,
                       mapping_proteins_to_proteoforms, proteoforms.str_to_int,
                       adjacency_list_proteoforms);
}

umsb convertSets(const umsb &traits_to_original_entities,
                 const vs &index_to_original_entities,
                 const ummss &mapping,
                 const umsi &result_entities_to_index,
                 const ummss &adjacency_list_result_entities) {
    umsb traits_to_result_entities;
    for (const auto &trait_entry : traits_to_original_entities) {  // For each trait entry
        uss candidates;
        for (int I = 0; I < traits_to_original_entities.size(); I++) {  // For each member in this trait
            if (trait_entry.second[I]) {
                auto range = mapping.equal_range(index_to_original_entities[I]);  // For each mapping entity
                for (auto it = range.first; it != range.second; it++) {
                    candidates.insert(it->second);
                }
            }
        }

        // Keep only those connected to any of the other gene set members in the reference network
        for (const auto &candidate : candidates) {
            auto range = adjacency_list_result_entities.equal_range(candidate);
            for (auto it = range.first; it != range.second; ++it) {
                if (candidates.find(it->second) != candidates.end()) {
                    // Set that the candidate is in the new set
                    traits_to_result_entities[trait_entry.first][result_entities_to_index.at(candidate)].set();
                    break;
                }
            }
        }
    }

    return traits_to_result_entities;
}

//uss getGeneStrings(const bitset<PHEGENI_GENES>& gene_set, const bimap_str_int& genes) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (gene_set.test(I)) {
//			result.insert(genes.int_to_str[I]);
//		}
//	}
//	return result;
//}
//
//uss getProteinStrings(const bitset<PHEGENI_PROTEINS>& protein_set, const bimap_str_int& proteins) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (protein_set.test(I)) {
//			result.insert(proteins.int_to_str[I]);
//		}
//	}
//	return result;
//}
//
//uss getProteoformStrings(const bitset<PHEGENI_PROTEOFORMS>& proteoform_set, const bimap_str_int& proteoforms) {
//	uss result;
//	for (int I = 0; I < PHEGENI_GENES; I++) {
//		if (proteoform_set.test(I)) {
//			result.insert(proteoforms.int_to_str[I]);
//		}
//	}
//	return result;
//}
//