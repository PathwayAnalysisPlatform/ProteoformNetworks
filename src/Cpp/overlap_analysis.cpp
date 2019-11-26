#include "overlap_analysis.hpp"

void report_overlaps_with_ptms();

using namespace std;

void doOverlapAnalysis(
        std::string_view path_file_modules,
        std::string_view path_file_genes,
        std::string_view path_file_proteins,
        std::string_view path_file_proteoforms,
        std::string_view path_file_gene_interactions,
        std::string_view path_file_protein_interactions,
        std::string_view path_file_proteoform_interactions,
        std::string_view path_file_mapping_proteins_to_genes,
        std::string_view path_file_mapping_proteins_to_proteoforms,
        std::string path_reports) {

    std::cout << "Reading genes, proteins and proteoforms.\n\n";
    const auto[genes, proteins, proteoforms] = get_entities(path_file_genes,
                                                            path_file_proteins,
                                                            path_file_proteoforms);

    // Read traits and genes in Phegeni file. The genes in the PhegenI file are ignored,
    // the universe set of genes used is the one of Reactome
    std::cout << "Reading gene modules.\n\n";
    const auto[traits, phegen_genes] = loadPheGenIGenesAndTraits(path_file_modules, genes);

    // Create or read module files at the three levels: all in one, and single module files.
    std::cout << "Creating protein and proteoform modules.\n\n";
    std::string path_modules = path_reports;
    path_modules += "modules/";
    auto all_modules = get_or_create_modules(
            path_modules,
            path_file_modules,
            path_file_gene_interactions,
            path_file_mapping_proteins_to_genes,
            path_file_protein_interactions,
            path_file_mapping_proteins_to_proteoforms,
            path_file_proteoform_interactions,
            genes, proteins, proteoforms,
            traits);

    for (const auto &level : levels)
        std::cerr << level << " level: " << all_modules.at(level).group_to_members.size() << " traits\n";

    // Check variation in module sizes
    std::cout << "Calculating module size variation.\n\n";
    report_module_size_variation(path_reports, all_modules, traits, MIN_MODULE_SIZE, MAX_MODULE_SIZE);

    // Calculate scores with Overlap Similarity
    std::cout << "Reading interaction networks...\n";
    const std::map<const std::string, const vusi> interactions = {
            {levels.at(0), loadInteractionNetwork(path_file_gene_interactions, genes, true)},
            {levels.at(1), loadInteractionNetwork(path_file_protein_interactions, proteins, true)},
            {levels.at(2), loadInteractionNetwork(path_file_proteoform_interactions, proteoforms, true)}
    };
    std::cout << "Calculating module pairs overlap.\n\n";
    report_pairs_overlap_data(path_reports, all_modules, interactions, traits, MIN_MODULE_SIZE, MAX_MODULE_SIZE);

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
const std::map<const std::string, const modules>
get_or_create_modules(const std::string &path_modules,
                      std::string_view path_file_phegeni,
                      std::string_view path_file_gene_interactions,
                      std::string_view path_file_mapping_proteins_to_genes,
                      std::string_view path_file_protein_interactions,
                      std::string_view path_file_mapping_proteins_to_proteoforms,
                      std::string_view path_file_proteoform_interactions,
                      const bimap_str_int &genes,
                      const bimap_str_int &proteins,
                      const bimap_str_int &proteoforms,
                      const bimap_str_int &traits
) {
    std::string suffix = "_modules.tsv";
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

    return {
            {
                    levels.at(0), gene_modules},
            {
                    levels.at(1), protein_modules},
            {
                    levels.at(2), proteoform_modules}
    };
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
    const bimap_str_int genes = createBimap(path_file_reactome_genes, true);
    std::cout << "Loading proteins..." << std::endl;
    const bimap_str_int proteins = createBimap(path_file_reactome_proteins, true);
    std::cout << "Loading proteoforms..." << std::endl;
    const bimap_str_int proteoforms = createBimap(path_file_reactome_proteoforms, false);

    std::cout << "Entities loaded.\n";
    return {genes, proteins, proteoforms};
}

/**
 * Creates three files with the overlap data for genes, proteins and proteoforms.
 *
 * The files have one row for each pair of modules that overlap in at least one entity.
 * The overlap data consists of overlap size (number of shared entities), overlap coefficient, interface size (number
 * of connections between nodes in different modules, or from an overlap node to a non overlap node)
 *
 * @param path_out Directory for the output files
 * @param all_modules Map with the modules for genes, proteins and proteoforms
 * @param traits Bimap with module names
 */
void report_pairs_overlap_data(const std::string &path_out,
                               const std::map<const std::string, const modules> &all_modules,
                               const std::map<const std::string, const vusi> &interactions,
                               const bimap_str_int &traits,
                               const int min_module_size, const int max_module_size) {

    std::cout << "Calculating overlap sizes at " << "genes" << " level. " << std::endl;
    auto overlap_sizes_proteoforms = getScores(all_modules.at(levels[2]).group_to_members, getOverlapSize,
                                   min_module_size, max_module_size);

    for (const auto &level : levels) {
        std::cout << "Calculating overlap sizes at " << level << " level. " << std::endl;
        auto overlap_sizes = getScores(all_modules.at(level).group_to_members, getOverlapSize,
                                       overlap_sizes_proteoforms);

        std::cout << "Calculating overlap similarity at " << level << " level." << std::endl;
        auto overlap_coefficients = getScores(all_modules.at(level).group_to_members, getOverlapSimilarity,
                                              overlap_sizes_proteoforms);

        std::cout << "Calculating node interface sizes at " << level << " level. " << std::endl;
        auto interface_sizes_nodes = getScores(all_modules.at(level).group_to_members, interactions.at(level),
                                         calculate_interface_size_nodes, overlap_sizes_proteoforms);

        std::cout << "Calculating edge interface sizes at " << level << " level. " << std::endl;
        auto interface_sizes_edges = getScores(all_modules.at(level).group_to_members, interactions.at(level),
                                               calculate_interface_size_nodes, overlap_sizes_proteoforms);

        std::cout << "Writing results to file." << std::endl;
        std::vector<std::string> features_labels = {"OVERLAP_COEFFICIENT", "NODE_OVERLAP_SIZE", "NODE_INTERFACE_SIZE", "EDGE_INTERFACE_SIZE" };
        std::vector<pair_map<double>> features = {overlap_coefficients, overlap_sizes, interface_sizes_nodes, interface_sizes_edges};
        std::string file_output = path_out + "pairs_overlap_data_" + level + ".tsv";

        writeScores(traits, all_modules.at(level), features_labels, features, file_output);
        std::cout << "Finished with " << level << "\n\n";
    }
}

void report_module_size_variation(std::string_view path_reports,
                                  const std::map<const std::string, const modules> &all_modules,
                                  const bimap_str_int &traits, const int min_module_size, const int max_module_size) {

    std::map<const std::string, um<int, int>> sizes = calculate_and_report_sizes(path_reports, all_modules, traits);

    // Check sizes of modules
    // Calculate how many empty modules exist, for genes, proteins and proteoforms
    std::cout << "Total number of possible modules is: " << traits.int_to_str.size() << std::endl;
    for (const auto &level : levels)
        std::cout << "Number of non-empty modules for " << level << ": " << sizes.at(level).size() << std::endl;

    std::ofstream output(path_reports.data() + static_cast<std::string>("sizes_variation.tsv"));

    if (!output.is_open()) {
        std::string message = "Cannot open report file at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    output << "TRAIT\tGENES_TO_PROTEINS\tPROTEINS_TO_PROTEOFORMS\n";
    for (const auto &module_entry : sizes.at(levels.at(0))) {
        output << module_entry.first << "\t";

        // Difference genes --> proteins
        if (hasKey(sizes.at(levels.at(1)), module_entry.first)) {
            output << sizes.at(levels.at(1)).at(module_entry.first) - sizes.at(levels.at(0)).at(module_entry.first);
        } else {
            output << -sizes.at(levels.at(0)).at(module_entry.first);
        }
        output << "\t";

        // Difference proteins --> proteoforms
        if (hasKey(sizes.at(levels.at(1)), module_entry.first)) {
            if (hasKey(sizes.at(levels.at(2)), module_entry.first)) {
                output << sizes.at(levels.at(2)).at(module_entry.first) - sizes.at(levels.at(1)).at(module_entry.first);
            } else {
                output << -sizes.at(levels.at(1)).at(module_entry.first);
            }
        } else {
            if (hasKey(sizes.at(levels.at(2)), module_entry.first)) {
                output << sizes.at(levels.at(2)).at(module_entry.first);
            } else {
                output << 0;
            }
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
    // For all the pairs of traits that overlap at the proteoform level

    // If overlapping proteoforms are modified

    // Report that pair of modules with the modifications of the proteins
}

// TODO: Test example of overlap pairs Adiponectin and Glomerular Filtration Rate. Check why the protein P15692 is not in the proteoform module of Glomerular Filtration Rate.
