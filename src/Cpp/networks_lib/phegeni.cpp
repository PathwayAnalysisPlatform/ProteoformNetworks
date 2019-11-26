#include "phegeni.hpp"

using namespace std;

// Read the genes and traits of PheGenI data file.
// Creates two bimaps: genes and traits
// The gene list are those genes in the dataset that are also contained in the acceptable gene list.
module_bimaps loadPheGenIGenesAndTraits(
        string_view path_file_phegeni,
        const bimap_str_int &acceptable_genes) {

    std::cout << "Loading PhegenI genes and traits." << std::endl;

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

    return {phegeni_traits, phegeni_genes};
}

// Read Phegeni trait modules with genes as members. Only gene members also in the acceptable gene list.
// Creates two mappings to represent the Trait modules with genes as members:
// - trait strings to gene bitsets.
// - gene strings to trait bitsets
// The two mappings define the trait modules (also called, phenotype modules or disease modules)
modules createAndLoadPheGenIGeneModules(std::string_view path_file_phegeni, const bimap_str_int &genes,
                                        const bimap_str_int &traits, const std::string &path_file_gene_interactions,
                                        const std::string &path_modules, const std::string &suffix,
                                        bool keep_disconnected_nodes) {
    ifstream file_phegen(path_file_phegeni.data());
    string line, field, trait, gene, gene2;
    string p_value_str;
    long double p_value;

    if (!file_phegen.is_open()) {
        std::string message = "Cannot open path_file_phegeni at ";
        std::string function = __FUNCTION__;
        throw runtime_error(message + function);
    }

    // Initialize modules
    modules modules;
    for (int trait_index = 0; trait_index < traits.int_to_str.size(); trait_index++)
        modules.group_to_members.push_back(base::dynamic_bitset<>(genes.int_to_str.size()));
    for (int gene_index = 0; gene_index < genes.int_to_str.size(); gene_index++)
        modules.member_to_groups.push_back(base::dynamic_bitset<>(traits.int_to_str.size()));

    std::cerr << "Initialized modules:\n"
              << "Genes to traits: " << modules.member_to_groups.size()
              << "\t Traits bitsets size: " << modules.member_to_groups.begin()->size()
              << "\nTraits to genes: " << modules.group_to_members.size()
              << "\t Genes bitsets size: " << modules.group_to_members.begin()->size() << "\n";

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

        if (hasKey(traits.str_to_int, trait) && hasKey(genes.str_to_int, gene)) {
            modules.group_to_members[traits.str_to_int.at(trait)][genes.str_to_int.at(gene)].set();
            modules.member_to_groups[genes.str_to_int.at(gene)][traits.str_to_int.at(trait)].set();
        }
        if (hasKey(traits.str_to_int, trait) && hasKey(genes.str_to_int, gene2)) {
            modules.group_to_members[traits.str_to_int.at(trait)][genes.str_to_int.at(gene2)].set();
            modules.member_to_groups[genes.str_to_int.at(gene2)][traits.str_to_int.at(trait)].set();
        }
    }

    // Load interaction network
    std::cout << "Reading gene interaction network..." << std::endl;

    if (!keep_disconnected_nodes) {
        std::cout << "Removing disconnected genes from modules..." << std::endl;
        auto interactions = loadInteractionNetwork(path_file_gene_interactions, genes, true);
        removeDisconnectedMembers(modules, traits, genes, interactions);
        writeModulesSingleFile(path_modules, "genes", suffix, modules, traits, genes);
        writeModulesManyFiles(path_modules, "genes", ".tsv", modules, traits, genes, interactions);
    }

    cerr << "Number of traits with gene members as bitset: " << modules.group_to_members.size() << "\n";
    cerr << "Number of genes with traits they belong as bitset: " << modules.member_to_groups.size() << "\n";

    return modules;
}

// Creates modules of proteins and proteoforms.
// Converts the modules from gene --> protein or protein --> proteoform, and then removes the disconnected vertices
// It applies the logic that, some nodes should not be in the network module, because there are no connection to others
modules createAndLoadPheGenIModules(const modules &prev_modules, const bimap_str_int &prev_entities,
                                    const bimap_str_int &entities, const bimap_str_int &traits, const ummss &mapping,
                                    const std::string &path_file_entity_interactions, const std::string &path_modules,
                                    const std::string &level, const std::string &suffix, bool keep_disconnected_nodes) {
    // Convert gene modules to protein modules using the mapping from entities to genes
    modules result_modules = convertModulesWithMapping(prev_modules, prev_entities, entities, traits, mapping);

    if (!keep_disconnected_nodes) {
        auto interactions = loadInteractionNetwork(path_file_entity_interactions, entities, true);
        removeDisconnectedMembers(result_modules, traits, entities, interactions);
        writeModulesSingleFile(path_modules, level, suffix, result_modules, traits, entities);
        writeModulesManyFiles(path_modules, level, ".tsv", result_modules, traits, entities, interactions);
    }

    return result_modules;
}

// Note that mapping uses strings, whereas modules use the indexes
// Therefore, the bimaps are necessary to jumpt between indexes and strings
modules convertModulesWithMapping(
        const modules &original_modules,
        const bimap_str_int &original_entities,
        const bimap_str_int &destination_entities,
        const bimap_str_int &traits,
        const ummss &mapping) {
    modules new_modules;

    for (int group_index = 0; group_index < traits.int_to_str.size(); group_index++)
        new_modules.group_to_members.push_back(base::dynamic_bitset<>(destination_entities.int_to_str.size()));

    for (int entity_index = 0; entity_index < destination_entities.int_to_str.size(); entity_index++)
        new_modules.member_to_groups.push_back(base::dynamic_bitset<>(traits.int_to_str.size()));

    for (int group_index = 0; group_index < traits.int_to_str.size(); group_index++) {
        for (int member_index = 0; member_index < original_entities.int_to_str.size(); member_index++) {
            if (original_modules.group_to_members[group_index][member_index]) {    // If it is really a member
                std::string original_entity = original_entities.int_to_str[member_index];

                auto range = mapping.equal_range(original_entity);
                // For each destination entity string of the original entity
                for (auto it = range.first; it != range.second; it++) {
                    // Check the mapped entity exists in the list of known acceptable destination entities
                    std::string destination_entity = it->second;
                    if (hasKey(destination_entities.str_to_int, destination_entity)) {
                        int destination_entity_index = destination_entities.str_to_int.at(destination_entity);

                        // Set the entities as members of the trait member set
                        new_modules.group_to_members[group_index][destination_entity_index].set();

                        // Set the traits as owners of the entity owner set
                        new_modules.member_to_groups[destination_entity_index][group_index].set();
                    } else {
//                        std::cerr << original_entity << " --> ?" << std::endl;
                    }
                }
            }
        }
    }

    return new_modules;
}

//

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