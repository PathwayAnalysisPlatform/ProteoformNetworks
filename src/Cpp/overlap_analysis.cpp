#include "overlap_analysis.hpp"

void report_overlaps_with_ptms();

using namespace std;

void doOverlapAnalysis(
        std::string_view path_file_phegeni,
        std::string_view path_file_reactome_genes,
        std::string_view path_file_reactome_proteins,
        std::string_view path_file_reactome_proteoforms,
        std::string_view path_file_mapping_proteins_to_genes,
        std::string_view path_file_mapping_proteins_to_proteoforms,
        std::string_view path_file_gene_interactions,
        std::string_view path_file_protein_interactions,
        std::string_view path_file_proteoform_interactions,
        std::string path_reports,
        std::string_view path_modules) {

    const auto[genes, proteins, proteoforms] = get_entities(path_file_reactome_genes,
                                                            path_file_reactome_proteins,
                                                            path_file_reactome_proteoforms);

    // Read traits and genes in Phegeni file. The genes in the PhegenI file are ignored,
    // the universe set of genes used is the one of Reactome
    const auto[traits, phegen_genes] = loadPheGenIGenesAndTraits(path_file_phegeni, genes);

    // Create or read module files at the three levels: all in one, and single module files.
    const auto[gene_modules, protein_modules, proteoform_modules] = get_or_create_modules(
            path_modules.data(),
            path_file_phegeni,
            path_file_gene_interactions,
            path_file_mapping_proteins_to_genes,
            path_file_protein_interactions,
            path_file_mapping_proteins_to_proteoforms,
            path_file_proteoform_interactions,
            genes, proteins, proteoforms,
            traits);

    std::cerr << "Gene level: " << gene_modules.group_to_members.size() << " traits\n";
    std::cerr << "Protein level: " << protein_modules.group_to_members.size() << " traits\n";
    std::cerr << "Proteoform level: " << proteoform_modules.group_to_members.size() << " traits\n";

    // Check variation in module sizes
    report_module_size_variation(path_reports, gene_modules, protein_modules, proteoform_modules, traits);

    // Calculate scores with Overlap Similarity
    const auto scores_overlap_similarity = get_scores(path_reports, getOverlapSimilarity, "overlap_similarity",
                                                      gene_modules, protein_modules, proteoform_modules, traits);

    report_node_overlap_reduction_examples(path_reports, "overlap_similarity", scores_overlap_similarity, traits);

    report_connecting_edges_variation_examples(path_reports, scores_overlap_similarity);

    report_overlap_only_ptms(path_reports, scores_overlap_similarity, traits);

    // Calculate scores with Jaccard index
//    const auto scores_jaccard_index = get_scores(path_scores, getJaccardSimilarity, "jaccard_index",
//                                                 gene_modules, protein_modules, proteoform_modules);

    // Check if there are modules overlapping at the proteoform level, only with modified proteins

    std::cout << "Finished with the Overlap analysis";
}

// Check if there are pairs of modules overlapping at one level but not in another
// Selects the pair of modules that had a higher overlap score at gene level than at protein or proteoform level
// Precondition: The modules for genes, proteins and proteoforms, must be the same, even when at some level there
// are empty modules. This allows comparing module indexes faster.
void report_node_overlap_reduction_examples(std::string path_scores, std::string label, const score_maps &scores,
                                            const bimap_str_int &traits) {

    std::ofstream output(path_scores.data() + label + "_score_variation_examples.tsv");

    if (!output.is_open()) {
        std::string message = "Cannot open report file at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    output << "TRAIT1\tTRAIT2\tSCORE_GENES\tSCORE_PROTEINS\tSCORE_PROTEOFORMS\t"
           << "GENES_TO_PROTEINS\tPROTEINS_TO_PROTEOFORMS\n";
    for (const auto &score_entry : scores.gene_scores) {
        // The score maps use the indexes of the traits, then
        if ((scores.protein_scores.at(score_entry.first) < scores.gene_scores.at(score_entry.first))
            || (scores.proteoform_scores.at(score_entry.first) < scores.gene_scores.at(score_entry.first))) {
            auto trait_index_pair = score_entry.first;
            output << traits.int_to_str[trait_index_pair.first] << '\t'
                   << traits.int_to_str[trait_index_pair.second] << "\t"
                   << scores.gene_scores.at(trait_index_pair) << "\t"
                   << scores.protein_scores.at(trait_index_pair) << "\t"
                   << scores.proteoform_scores.at(trait_index_pair) << "\t"
                   << scores.protein_scores.at(trait_index_pair) - scores.gene_scores.at(trait_index_pair) << "\t"
                   << scores.proteoform_scores.at(trait_index_pair) - scores.protein_scores.at(trait_index_pair)
                   << "\n";
        }
    }

    output.close();
}

// Create or read module files at the three levels: all in one, and single module files.
// There should be three files with all modules, one for genes, one for proteins and one for protoeforms.
// There should be three more files for each single trait in an individual file.
get_modules_result get_or_create_modules(std::string path_modules,
                                         std::string_view path_file_phegeni,
                                         std::string_view path_file_gene_interactions,
                                         std::string_view path_file_mapping_proteins_to_genes,
                                         std::string_view path_file_protein_interactions,
                                         std::string_view path_file_mapping_proteins_to_proteoforms,
                                         std::string_view path_file_proteoform_interactions,
                                         const bimap_str_int &genes,
                                         const bimap_str_int &proteins,
                                         const bimap_str_int &proteoforms,
                                         const bimap_str_int &traits) {
    std::string suffix = ".tsv";
    modules gene_modules, protein_modules, proteoform_modules;

    std::string all_traits_genes_file_name = path_modules + "genes" + suffix;
    std::string all_traits_proteins_file_name = path_modules + "proteins" + suffix;
    std::string all_traits_proteoforms_file_name = path_modules + "proteoforms" + suffix;

    if (file_exists(all_traits_genes_file_name)) {
        std::cout << "Reading gene modules." << std::endl;
        gene_modules = loadModules(all_traits_genes_file_name, traits, genes).entity_modules;
    } else {
        std::cout << "Creating gene modules." << std::endl;
        gene_modules = createAndLoadPheGenIGeneModules(path_file_phegeni, genes, traits,
                                                       path_file_gene_interactions.data(),
                                                       path_modules, suffix);
    }

    if (file_exists(all_traits_proteins_file_name)) {
        std::cout << "Reading protein modules." << std::endl;
        protein_modules = loadModules(all_traits_proteins_file_name, traits, proteins).entity_modules;
    } else {
        std::cout << "Creating protein modules." << std::endl;
        bidirectional_mapping mapping_genes_to_proteins = readMapping(path_file_mapping_proteins_to_genes,
                                                                      true, false);
        protein_modules = createAndLoadPheGenIModules(gene_modules, genes, proteins, traits,
                                                      mapping_genes_to_proteins.second_to_first,
                                                      path_file_protein_interactions.data(),
                                                      path_modules, "proteins", suffix);
    }

    if (file_exists(all_traits_proteoforms_file_name)) {
        std::cout << "Reading proteoform modules." << std::endl;
        proteoform_modules = loadModules(all_traits_proteoforms_file_name, traits, proteoforms).entity_modules;
    } else {
        std::cout << "Creating proteoform modules." << std::endl;
        bidirectional_mapping mapping_proteins_to_proteoforms = readMapping(path_file_mapping_proteins_to_proteoforms,
                                                                            true, true);
        proteoform_modules = createAndLoadPheGenIModules(protein_modules, proteins, proteoforms, traits,
                                                         mapping_proteins_to_proteoforms.second_to_first,
                                                         path_file_proteoform_interactions.data(),
                                                         path_modules, "proteoforms", suffix);
    }

    return {gene_modules, protein_modules, proteoform_modules};
}


void printMembers(std::ostream &output, const uss &members) {
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

// Read Reactome genes, proteins and proteoforms. Take them as universe set of acceptable identifiers.
get_entities_result get_entities(std::string_view path_file_reactome_genes,
                                 std::string_view path_file_reactome_proteins,
                                 std::string_view path_file_reactome_proteoforms) {

    std::cout << "Loading genes..." << std::endl;
    const bimap_str_int genes = createBimap(path_file_reactome_genes);
    std::cout << "Loading proteins..." << std::endl;
    const bimap_str_int proteins = createBimap(path_file_reactome_proteins);
    std::cout << "Loading proteoforms..." << std::endl;
    const bimap_str_int proteoforms = createBimap(path_file_reactome_proteoforms, true, 0, 6);

    return {genes, proteins, proteoforms};
}

score_maps get_scores(std::string path_scores,
                      std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>)> scoring,
                      std::string label,
                      const modules &gene_modules,
                      const modules &protein_modules,
                      const modules &proteoform_modules,
                      const bimap_str_int &traits) {

    std::cout << "Calculating overlap similarity for Genes." << std::endl;
    auto gene_scores = getScores(gene_modules.group_to_members, scoring);
    writeScores(traits, gene_modules, gene_scores, path_scores + "scores_genes_" + label + ".tsv");

    std::cout << "Calculating overlap similarity for Proteins." << std::endl;
    auto protein_scores = getScores(protein_modules.group_to_members, scoring);
    writeScores(traits, protein_modules, protein_scores, path_scores + "scores_proteins_" + label + ".tsv");

    std::cout << "Calculating overlap similarity for Proteoforms." << std::endl;
    auto proteoform_scores = getScores(proteoform_modules.group_to_members, scoring);
    writeScores(traits, proteoform_modules, proteoform_scores, path_scores + "scores_proteoforms_" + label + ".tsv");

    return {gene_scores, protein_scores, proteoform_scores};
}

void report_module_size_variation(std::string_view path_reports, const modules &gene_modules,
                                  const modules &protein_modules, const modules &proteoform_modules,
                                  const bimap_str_int &traits) {
    const auto sizes_genes = calculate_and_report_sizes(path_reports, "genes", gene_modules, traits);
    const auto sizes_proteins = calculate_and_report_sizes(path_reports, "proteins", protein_modules, traits);
    const auto sizes_proteoforms = calculate_and_report_sizes(path_reports, "proteoforms", proteoform_modules, traits);

    // Check sizes of modules
    // Calculate how many empty modules exist, for genes, proteins and proteoforms
    std::cout << "Total number of possible modules is: " << traits.int_to_str.size() << std::endl;
    std::cout << "Number of non-empty gene modules: " << sizes_genes.size() << std::endl;
    std::cout << "Number of non-empty protein modules: " << sizes_proteins.size() << std::endl;
    std::cout << "Number of non-empty proteoform modules: " << sizes_proteoforms.size() << std::endl;

    std::ofstream output(path_reports.data() + static_cast<std::string>("sizes_variation.tsv"));

    if (!output.is_open()) {
        std::string message = "Cannot open report file at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    output << "TRAIT\tGENES_TO_PROTEINS\tPROTEINS_TO_PROTEOFORMS\tGENES_TO_PROTEOFORMS\n";
    for (const auto &module_entry : sizes_genes) {
        output << module_entry.first << "\t";

        // Difference genes --> proteins
        if (hasKey(sizes_proteins, module_entry.first)) {
            output << sizes_proteins.at(module_entry.first) - sizes_genes.at(module_entry.first);
        } else {
            output << -sizes_genes.at(module_entry.first);
        }
        output << "\t";

        // Difference proteins --> proteoforms
        if (hasKey(sizes_proteins, module_entry.first)) {
            if (hasKey(sizes_proteoforms, module_entry.first)) {
                output << sizes_proteoforms.at(module_entry.first) - sizes_proteins.at(module_entry.first);
            } else {
                output << -sizes_proteins.at(module_entry.first);
            }
        } else {
            if (hasKey(sizes_proteoforms, module_entry.first)) {
                output << sizes_proteoforms.at(module_entry.first);
            } else {
                output << 0;
            }
        }
        output << "\t";

        // Difference genes --> proteoforms
        if (hasKey(sizes_proteoforms, module_entry.first)) {
            output << sizes_proteoforms.at(module_entry.first) - sizes_genes.at(module_entry.first);
        } else {
            output << -sizes_genes.at(module_entry.first);
        }
        output << "\n";
    }
    output.close();
}

void report_connecting_edges_variation_examples(std::string path_reports, const score_maps scores) {

    std::ofstream output(path_reports.data() + static_cast<std::string>("bridges_variation.tsv"));

    if (!output.is_open()) {
        std::string message = "Cannot open report file at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    // Check those pairs of modules which overlap at the gene level
    for (const auto &pair_entry : scores.gene_scores) {
        std::pair<int, int> pair_names = pair_entry.first;
        double pair_score = pair_entry.second;
        if (pair_score > 0) {
            // Given two modules m1 and m2, and vertices v1 in m1 only, and v2 in m2 only.
            // A bridge edge is an edge connecting v1 and v2. An edge connecting the two modules.
            // Count the number of bridge edges at each level


            // Calculate the ratio of change between each level
        }
    }


}

void report_overlap_only_ptms(std::string string, const score_maps maps, const bimap_str_int anInt) {
    //TODO: Implement this function
}

