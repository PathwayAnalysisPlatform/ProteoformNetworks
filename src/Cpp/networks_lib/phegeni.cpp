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
trait_modules loadPheGenIGeneModules(
        string_view path_file_phegeni,
        string_view path_file_genes) {
    ifstream file_phegen(path_file_phegeni.data());
    string line, field, trait, gene, gene2;
    string p_value_str;
    long double p_value;

    if (!file_phegen.is_open()) {
        std::string message = "Cannot open path_file_phegeni at ";
        std::string function = __FUNCTION__;
        throw runtime_error(message + function);
    }

    // Read Reactome genes. Take them as all acceptable gene names.
    const bimap_str_int acceptable_genes = createBimap(path_file_genes); // Gene names --> str_to_int

    // Read traits and genes in Phegeni file
    const auto[phegeni_genes, phegeni_traits] = loadPheGenIGenesAndTraits(path_file_phegeni, acceptable_genes);

    // Initialize modules
    trait_modules modules;
    for (const auto &gene : phegeni_genes.int_to_str) {
        modules.entities_to_traits[gene] = base::dynamic_bitset<>(phegeni_traits.int_to_str.size());
    }
    for (const auto &trait : phegeni_traits.int_to_str) {
        modules.traits_to_entities[trait] = base::dynamic_bitset<>(phegeni_genes.int_to_str.size());
    }

    std::cout << "Initialized modules:\n"
              << "Genes to traits: " << modules.entities_to_traits.size()
              << "\t Traits bitsets size: " << modules.entities_to_traits.begin()->second.size()
              << "\nTraits to genes: " << modules.traits_to_entities.size()
              << "\t Genes bitsets size: " << modules.traits_to_entities.begin()->second.size() << "\n";

    // Read members of each module
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
        getline(file_phegen,
                line);               // Skip header line leftoever: Source,	PubMed,	Analysis ID,	Study ID,	Study Name

        if (acceptable_genes.str_to_int.find(gene) != acceptable_genes.str_to_int.end()) {
            modules.traits_to_entities[trait][phegeni_genes.str_to_int.at(gene)].set();
            modules.entities_to_traits[gene][phegeni_traits.str_to_int.at(trait)].set();
        }
        if (acceptable_genes.str_to_int.find(gene2) != acceptable_genes.str_to_int.end()) {
            modules.traits_to_entities[trait][phegeni_genes.str_to_int.at(gene2)].set();
            modules.entities_to_traits[gene2][phegeni_traits.str_to_int.at(trait)].set();
        }
    }

    cerr << "Number of traits with gene members as bitset: " << modules.traits_to_entities.size() << "\n";
    cerr << "Number of genes with traits they belong as bitset: " << modules.entities_to_traits.size() << "\n";

    return modules;
}

trait_modules createPheGenIProteinModules(const trait_modules &gene_modules,
                                          const bimap_str_int &genes,
                                          const bimap_str_int &proteins,
                                          std::string_view path_file_proteins_to_genes,
                                          std::string_view path_file_protein_edges) {
    // Convert gene modules to protein modules using the mapping from proteins to genes
    entity_mapping mapping = readMapping(path_file_proteins_to_genes);
    trait_modules protein_modules = convertModulesWithMapping(gene_modules, genes, proteins, mapping);

    // Load interaction network
    ummss protein_interactions;
//    const auto[adjacency_list_proteins, adjacency_list_proteoforms] = loadReactomeNetworks(path_file_protein_edges,
//                                                                                           path_file_proteoform_edges);

    // For each module, discard the members which are not connected to other members in the interaction network

    return protein_modules;
}

trait_modules convertModulesWithMapping(
        const trait_modules &original_modules,
        const bimap_str_int &original_entities,
        const bimap_str_int &destination_entities,
        const entity_mapping &mapping) {
    trait_modules modules;
//    for (const auto &trait_entry : traits_to_original_entities) {  // For each trait entry
//        uss candidates;
//        for (int I = 0; I < traits_to_original_entities.size(); I++) {  // For each member in this trait
//            if (trait_entry.second[I]) {
//                auto range = mapping.equal_range(index_to_original_entities[I]);  // For each mapping entity
//                for (auto it = range.first; it != range.second; it++) {
//                    candidates.insert(it->second);
//                }
//            }
//        }
//
//        // Keep only those connected to any of the other gene set members in the reference network
//        for (const auto &candidate : candidates) {
//            auto range = adjacency_list_result_entities.equal_range(candidate);
//            for (auto it = range.first; it != range.second; ++it) {
//                if (candidates.find(it->second) != candidates.end()) {
//                    // Set that the candidate is in the new set
//                    traits_to_result_entities[trait_entry.first][result_entities_to_index.at(candidate)].set();
//                    break;
//                }
//            }
//        }
//    }

    return modules;
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