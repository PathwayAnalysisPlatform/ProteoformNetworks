#include "overlap_analysis.hpp"

void get_modules();

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
        std::string path_scores,
        std::string_view path_modules) {

    const auto[genes, proteins, proteoforms] = get_entities(path_file_reactome_genes,
                                                            path_file_reactome_proteins,
                                                            path_file_reactome_proteoforms);

    // Read traits and genes in Phegeni file. The genes in the PhegenI file are ignored,
    // the universe set of genes used is the one of Reactome
    const auto[traits, phegen_genes] = loadPheGenIGenesAndTraits(path_file_phegeni, genes);

    const auto[gene_modules, protein_modules, proteoform_modules] = get_or_create_modules(
            path_modules.data(), path_file_phegeni, path_file_gene_interactions,
            path_file_mapping_proteins_to_genes, path_file_protein_interactions,
            path_file_mapping_proteins_to_proteoforms, path_file_proteoform_interactions,
            genes, proteins, proteoforms, traits);

    // Calculate scores with Overlap Similarity
    const auto scores_overlap_similarity = get_scores(path_scores, getOverlapSimilarity, "overlap_similarity",
                                                      gene_modules, protein_modules, proteoform_modules);

    // Calculate scores with Jaccard index
    const auto scores_jaccard_index = get_scores(path_scores, getJaccardSimilarity, "jaccard_index",
                                                 gene_modules, protein_modules, proteoform_modules);

    // Calculate how many empty modules exist, for genes, proteins and proteoforms

    // Check if the overlap changed from genes to proteins to proteoforms, for non-empty module pairs

    // Calculate the frequency of score values for non-empty module pairs

    // Check if there are modules overlapping at one level but not in another

    // Check if there are modules overlapping at the proteoform level, only with modified proteins

    std::cout << "Finished with the Overlap analysis";
}

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
    std::string suffix = "_modules.tsv";
    modules gene_modules, protein_modules, proteoform_modules;

    if (file_exists(path_modules + "gene" + suffix)) {
        std::cout << "Reading gene modules." << std::endl;
        gene_modules = loadModules(path_modules + "gene" + suffix).entity_modules;
    } else {
        std::cout << "Creating gene modules." << std::endl;
        gene_modules = loadPheGenIGeneModules(path_file_phegeni, genes, traits, path_file_gene_interactions);
        writeModules(path_modules + "gene" + suffix, gene_modules, traits, genes);
    }

    if (file_exists(path_modules + "protein" + suffix)) {
        std::cout << "Reading protein modules." << std::endl;
        protein_modules = loadModules(path_modules + "protein" + suffix).entity_modules;
    } else {
        std::cout << "Creating protein modules." << std::endl;
        bidirectional_mapping mapping_genes_to_proteins = readMapping(path_file_mapping_proteins_to_genes, true,
                                                                      false);
        protein_modules = createPheGenIModules(gene_modules, genes, proteins, traits,
                                               mapping_genes_to_proteins.second_to_first,
                                               path_file_protein_interactions);
        writeModules(path_modules + "protein" + suffix, protein_modules, traits, proteins);
    }

    if (file_exists(path_modules + "proteoform" + suffix)) {
        std::cout << "Reading proteoform modules." << std::endl;
        proteoform_modules = loadModules(path_modules + "proteoform" + suffix).entity_modules;
    } else {
        std::cout << "Creating proteoform modules." << std::endl;
        bidirectional_mapping mapping_proteins_to_proteoforms = readMapping(path_file_mapping_proteins_to_proteoforms,
                                                                            true,
                                                                            true);
        proteoform_modules = createPheGenIModules(protein_modules, proteins, proteoforms, traits,
                                                  mapping_proteins_to_proteoforms.second_to_first,
                                                  path_file_proteoform_interactions);

        writeModules(path_modules + "proteoform" + suffix, proteoform_modules, traits, proteoforms);
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

get_scores_result get_scores(std::string path_scores,
                             std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>)> scoring,
                             std::string label,
                             const modules &gene_modules,
                             const modules &protein_modules,
                             const modules &proteoform_modules) {

    std::cout << "Calculating overlap similarity for Genes." << std::endl;
    auto gene_scores = getScores(gene_modules.group_to_members, scoring);
    writeScores(gene_modules.group_to_members, gene_scores,
                path_scores + "scores_genes_" + label + ".tsv");

    std::cout << "Calculating overlap similarity for Proteins." << std::endl;
    auto protein_scores = getScores(protein_modules.group_to_members, scoring);
    writeScores(protein_modules.group_to_members, protein_scores,
                path_scores + "scores_proteins_" + label + ".tsv");

    std::cout << "Calculating overlap similarity for Proteoforms." << std::endl;
    auto proteoform_scores = getScores(proteoform_modules.group_to_members, scoring);
    writeScores(proteoform_modules.group_to_members, proteoform_scores,
                path_scores + "scores_proteoforms_" + label + ".tsv");

    return {gene_scores, protein_scores, proteoform_scores};
}

void calculate_statistics() {
    // Calculate how many empty modules exist, for genes, proteins and proteoforms

    // Check if the overlap changed from genes to proteins to proteoforms, for non-empty module pairs

    // Calculate the frequency of score values for non-empty module pairs

    // Check if there are modules overlapping at one level but not in another

    // Check if there are modules overlapping at the proteoform level, only with modified proteins
}

