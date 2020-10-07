#include "overlap_analysis.hpp"

void report_overlaps_with_ptms();

using namespace std;

std::string GetExeFileName() {
    char buffer[MAX_PATH];
    GetModuleFileName(NULL, buffer, MAX_PATH);
    return std::string(buffer);
}

std::string GetExePath() {
    std::string f = GetExeFileName();
    return f.substr(0, f.find_last_of("\\/"));
}

void doOverlapAnalysis(
        std::string_view path_file_modules,
        std::string_view path_file_genes,
        std::string_view path_file_proteins,
        std::string_view path_file_proteoforms,
        std::string_view path_file_small_molecules,
        std::string_view path_file_gene_interactions,
        std::string_view path_file_protein_interactions,
        std::string_view path_file_proteoform_interactions,
        std::string_view path_file_mapping_proteins_to_genes,
        std::string_view path_file_mapping_proteins_to_proteoforms,
        std::string path_reports) {

//    std::cout << "current directory is " << GetExePath() << "\n";
//    std::cout << "executing file is " << GetExeFileName() << "\n";

    std::cout << "Reading genes, proteins and proteoforms.\n\n";
    const auto[genes, proteins, proteoforms, small_molecules] = get_entities(
            path_file_genes,
            path_file_proteins,
            path_file_proteoforms,
            path_file_small_molecules);


    // Read traits and genes in Phegeni file. The genes in the PhegenI file are ignored,
    // the universe set of genes used is the one of Reactome
    std::cout << "Reading gene All_modules.\n\n";
    const auto[traits, phegen_genes] = loadPheGenIGenesAndTraits(path_file_modules, genes);

    // Create or read module files at the three levels: all in one, and single module files.
    std::cout << "Creating All_modules\n\n";
    std::string path_modules = path_reports;

    auto all_modules = get_or_create_modules(
            path_file_modules,
            genes,
            proteins,
            proteoforms,
            small_molecules,
            path_file_gene_interactions,
            path_file_protein_interactions,
            path_file_proteoform_interactions,
            path_file_mapping_proteins_to_genes,
            path_file_mapping_proteins_to_proteoforms,
            path_reports,
            true);

//    int num_modules = all_modules.begin()->second.group_to_members.size();
//    for (const auto &entry : LEVELS) {
//        const char* level = entry.str;
//        std::cerr << "Num " << level << " All_modules: " << all_modules.at(level).group_to_members.size() << std::endl;
//        assert(num_modules == all_modules.at(level).group_to_members.size());
//    }
//    calculate_and_report_sizes(path_modules, all_modules, traits);
//
//    std::cout << "Creating All_modules removing disconnected nodes...\n\n";
//    path_modules = path_reports += "All_modules/";
//    all_modules = get_or_create_modules(
//            path_modules,
//            path_file_modules,
//            path_file_gene_interactions,
//            path_file_mapping_proteins_to_genes,
//            path_file_protein_interactions,
//            path_file_mapping_proteins_to_proteoforms,
//            path_file_proteoform_interactions,
//            genes, proteins, proteoforms, small_molecules,
//            traits, false);
//
//    num_modules = all_modules.begin()->second.group_to_members.size();
//    for (const auto &level : LEVELS) {
//        std::cerr << "Num " << level << " All_modules: " << all_modules.at(level).group_to_members.size() << std::endl;
//        assert(num_modules == all_modules.at(level).group_to_members.size());
//    }
//    calculate_and_report_sizes(path_reports, all_modules, traits);

    // Calculate scores with Overlap Similarity
//    std::cout << "Reading interaction networks...\n";
//    const std::map<const std::string, const vusi> interactions = {
//            {LEVELS.at(0), loadInteractionNetwork(path_file_gene_interactions, genes, true)},
//            {LEVELS.at(1), loadInteractionNetwork(path_file_protein_interactions, proteins, true)},
//            {LEVELS.at(2), loadInteractionNetwork(path_file_proteoform_interactions, proteoforms, true)}
//    };
//    std::cout << "Calculating module pairs overlap.\n\n";
//    report_pairs_overlap_data(path_reports, all_modules, interactions, traits, MIN_MODULE_SIZE, MAX_MODULE_SIZE);
//
//    std::cout << "Finished with the Overlap analysis";
}

// Check if there are pairs of All_modules overlapping at one level but not in another
// Selects the pair of All_modules that had a higher overlap score at gene level than at protein or proteoform level
// Precondition: The All_modules for genes, proteins and proteoforms, must be the same, even when at some level there
// are empty All_modules. This allows comparing module indexes faster.
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
// There should be three files with all All_modules, one for genes, one for proteins and one for proteoforms.
// There should be three more files for each single trait in an individual file.
std::map<const char *, const All_modules>
get_or_create_modules(std::string_view path_file_phegeni,
                      const entities_bimap &genes,
                      const entities_bimap &proteins,
                      const entities_bimap &proteoforms,
                      const entities_bimap &small_molecules,
                      std::string_view path_file_gene_interactions,
                      std::string_view path_file_protein_interactions,
                      std::string_view path_file_proteoform_interactions,
                      std::string_view path_file_mapping_proteins_to_genes,
                      std::string_view path_file_mapping_proteins_to_proteoforms,
                      const std::string &path_output,
                      bool keep_disconnected_nodes) {
    All_modules gene_modules, protein_modules, proteoform_modules;

    // Read traits and genes in Phegeni file. The returned genes in the PhegenI file are ignored,
    // the universe set of genes used is the one of Reactome
    std::cout << "Reading gene All_modules.\n\n";
    const auto[traits, phegen_genes] = loadPheGenIGenesAndTraits(path_file_phegeni, genes);

    gene_modules = get_or_create_gene_modules(path_file_phegeni,
                                              traits, genes, small_molecules, path_file_gene_interactions,
                                              path_output, keep_disconnected_nodes);

    // Proteins
//    modules_file_name = path_output + LEVELS.at(protein) + suffix;
//    if (file_exists(modules_file_name)) {
//        std::cout << "Reading protein All_modules." << std::endl;
//        protein_modules = loadModules(modules_file_name, traits, allEntities.proteins).entity_modules;
//    } else {
//        std::cout << "Creating protein All_modules." << std::endl;
//        bidirectional_mapping mapping_genes_to_proteins = readMapping(path_file_mapping_proteins_to_genes,
//                                                                      true, false);
//        protein_modules = createProteinOrProteoformModules(gene_modules, allEntities.genes, allEntities.proteins, traits,
//                                                           mapping_genes_to_proteins.second_to_first,
//                                                           path_file_protein_interactions.data(),
//                                                           path_output, LEVELS.at(protein), suffix,
//                                                           keep_disconnected_nodes);
//    }
//
//    // Proteoforms
//    modules_file_name = path_output + LEVELS.at(proteoform) + suffix;
//    if (file_exists(modules_file_name)) {
//        std::cout << "Reading proteoform All_modules." << std::endl;
//        proteoform_modules = loadModules(modules_file_name, traits, proteoforms).entity_modules;
//    } else {
//        std::cout << "Creating proteoform All_modules." << std::endl;
//        bidirectional_mapping mapping_proteins_to_proteoforms = readMapping(path_file_mapping_proteins_to_proteoforms,
//                                                                            true, true);
//        proteoform_modules = createProteinOrProteoformModules(protein_modules, proteins, proteoforms, traits,
//                                                              mapping_proteins_to_proteoforms.second_to_first,
//                                                              path_file_proteoform_interactions.data(),
//                                                              path_output, LEVELS.at(proteoform), suffix,
//                                                              keep_disconnected_nodes);
//    }

    return {{config::LEVELS.at(config::gene),       gene_modules},
            {config::LEVELS.at(config::protein),    protein_modules},
            {config::LEVELS.at(config::proteoform), proteoform_modules}};
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
All_entities get_entities(std::string_view file_genes,
                          std::string_view file_proteins,
                          std::string_view file_proteoforms,
                          std::string_view file_small_molecules) {

    std::cout << "Loading genes..." << std::endl;
    const bimap_str_int genes = createBimap(file_genes, false);
    std::cout << "Loading proteins..." << std::endl;
    const bimap_str_int proteins = createBimap(file_proteins, false);
    std::cout << "Loading proteoforms..." << std::endl;
    const bimap_str_int proteoforms = createBimap(file_proteoforms, false);
    std::cout << "Loading small molecules..." << std::endl;
    const bimap_str_int small_molecules = createBimap(file_small_molecules, false);

    std::cout << "Entities loaded.\n";
    return {genes, proteins, proteoforms, small_molecules};
}

/**
 * Creates three files with the overlap data for genes, proteins and proteoforms.
 *
 * The files have one row for each pair of All_modules that overlap in at least one entity.
 * The overlap data consists of overlap size (number of shared entities), overlap coefficient, interface size (number
 * of connections between nodes in different All_modules, or from an overlap node to a non overlap node)
 *
 * @param path_out Directory for the output files
 * @param all_modules Map with the All_modules for genes, proteins and proteoforms
 * @param traits Bimap with module names
 */
void report_pairs_overlap_data(const std::string &path_out,
                               const std::map<const std::string, const All_modules> &all_modules,
                               const std::map<const std::string, const vusi> &interactions,
                               const bimap_str_int &traits,
                               const int min_module_size, const int max_module_size) {

    std::cout << "Calculating overlap sizes at " << "genes" << " level. " << std::endl;
    auto overlap_sizes_proteoforms = getScores(all_modules.at(config::LEVELS.at(config::proteoform)).group_to_members, getOverlapSize,
                                               min_module_size, max_module_size);

    for (const auto &level : config::LEVELS) {

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
        std::vector<std::string> features_labels = {"OVERLAP_COEFFICIENT", "NODE_OVERLAP_SIZE", "NODE_INTERFACE_SIZE",
                                                    "EDGE_INTERFACE_SIZE"};
        std::vector<pair_map<double>> features = {overlap_coefficients, overlap_sizes, interface_sizes_nodes,
                                                  interface_sizes_edges};
        std::string file_output = path_out + "pairs_overlap_data_" + level + ".tsv";

        writeScores(traits, all_modules.at(level), features_labels, features, file_output);
        std::cout << "Finished with " << level << "\n\n";
    }
}

void report_module_size_variation(std::string_view path_reports,
                                  const std::map<const std::string, const All_modules> &all_modules,
                                  const bimap_str_int &traits, std::map<const std::string, um<int, int>> &sizes) {

    // Check sizes of All_modules
    // Calculate how many empty All_modules exist, for genes, proteins and proteoforms
    std::cout << "Total number of possible All_modules is: " << traits.int_to_str.size() << std::endl;
    std::cout << "Number of non-empty All_modules for " << config::LEVELS.at(config::gene) << ": " << sizes.at(config::LEVELS.at(config::gene)).size()
              << std::endl;
    std::cout << "Number of non-empty All_modules for " << config::LEVELS.at(config::protein) << ": " << sizes.at(config::LEVELS.at(config::protein)).size()
              << std::endl;
    std::cout << "Number of non-empty All_modules for " << config::LEVELS.at(config::proteoform) << ": "
              << sizes.at(config::LEVELS.at(config::proteoform)).size() << std::endl;

    std::ofstream output(path_reports.data() + static_cast<std::string>("sizes_variation.tsv"));

    if (!output.is_open()) {
        std::string message = "Cannot open report file at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    output << "TRAIT\tGENES_TO_PROTEINS\tPROTEINS_TO_PROTEOFORMS\n";
    for (const auto &module_entry : sizes.at(config::LEVELS.at(config::gene))) {
        output << module_entry.first << "\t";

        // Difference genes --> proteins
        if (hasKey(sizes.at(config::LEVELS.at(config::protein)), module_entry.first)) {
            output << sizes.at(config::LEVELS.at(config::protein)).at(module_entry.first) -
                      sizes.at(config::LEVELS.at(config::gene)).at(module_entry.first);
        } else {
            output << -sizes.at(config::LEVELS.at(config::gene)).at(module_entry.first);
        }
        output << "\t";

        // Difference proteins --> proteoforms
        if (hasKey(sizes.at(config::LEVELS.at(config::protein)), module_entry.first)) {
            if (hasKey(sizes.at(config::LEVELS.at(config::proteoform)), module_entry.first)) {
                output << sizes.at(config::LEVELS.at(config::proteoform)).at(module_entry.first) -
                          sizes.at(config::LEVELS.at(config::protein)).at(module_entry.first);
            } else {
                output << -sizes.at(config::LEVELS.at(config::protein)).at(module_entry.first);
            }
        } else {
            if (hasKey(sizes.at(config::LEVELS.at(config::proteoform)), module_entry.first)) {
                output << sizes.at(config::LEVELS.at(config::proteoform)).at(module_entry.first);
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

    // Check those pairs of All_modules which overlap at the gene level
    for (const auto &pair_entry : scores.gene_scores) {
        std::pair<int, int> pair_names = pair_entry.first;
        double pair_score = pair_entry.second;
        if (pair_score > 0) {
            // Given two All_modules m1 and m2, and vertices v1 in m1 only, and v2 in m2 only.
            // A bridge edge is an edge connecting v1 and v2. An edge connecting the two All_modules.
            // Count the number of bridge edges at each level


            // Calculate the ratio of change between each level
        }
    }
}

void report_overlap_only_ptms(std::string string, const score_maps maps, const bimap_str_int anInt) {
    //TODO: Implement this function
    // For all the pairs of traits that overlap at the proteoform level

    // If overlapping proteoforms are modified

    // Report that pair of All_modules with the modifications of the proteins
}

// TODO: Test example of overlap pairs Adiponectin and Glomerular Filtration Rate.
//  Check why the protein P15692 is not in the proteoform module of Glomerular Filtration Rate.

/*Performs the Overlap analysis calculations and creates files with the results.
 *
 * INPUT:
 * - [1] PheGenI gene sets file
 * - [2] Gene list (has header row)
 * - [3] Proteins list (has header row)
 * - [4] Proteoform list (does NOT have header row)
 * - [5] Gene interactions
 * - [6] Protein interactions
 * - [7] Proteoform interactions
 * - [8] Mapping from proteins to genes
 * - [9] Mapping from proteins to proteoforms
 * - [10] Output path
 *
 * */
int main(int argc, char *argv[]) try {
//	auto out = freopen("out.txt", "w", stdout);
//	auto err = freopen("err.txt", "w", stderr);

    if (argc < 11) {
        std::cerr << "Missing arguments. Expected: 10 arguments:\n\n"
                  << " * - [1] Modules file\n"
                  << " * - [2] Gene list\n"
                  << " * - [3] Proteins list\n"
                  << " * - [4] Proteoform list\n"
                  << " * - [5] Small molecules list\n"
                  << " * - [6] Gene interactions\n"
                  << " * - [7] Protein interactions\n"
                  << " * - [8] Proteoform interactions\n"
                  << " * - [9] Mapping from proteins to genes\n"
                  << " * - [10] Mapping from proteins to proteoforms\n"
                  << " * - [11] Output path";
        throw std::runtime_error("Missing arguments.");
        return 0;
    }

    doOverlapAnalysis(argv[1],
                      argv[2],
                      argv[3],
                      argv[4],
                      argv[5],
                      argv[6],
                      argv[7],
                      argv[8],
                      argv[9],
                      argv[10],
                      argv[11]);

}
catch (const std::exception &ex) {
    std::cout << ex.what() << "\n";
}