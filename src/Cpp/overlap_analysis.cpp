#include "overlap_analysis.hpp"

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

    // Read Reactome genes, proteins and proteoforms. Take them as set of acceptable identifiers.
    std::cout << "Loading genes..." << std::endl;
    const bimap_str_int genes = createBimap(path_file_reactome_genes);
    std::cout << "Loading proteins..." << std::endl;
    const bimap_str_int proteins = createBimap(path_file_reactome_proteins);
    std::cout << "Loading proteoforms..." << std::endl;
    const bimap_str_int proteoforms = createBimap(path_file_reactome_proteoforms, true, 0, 6);

    // Read traits and genes in Phegeni file
    std::cout << "Loading PhegenI genes and traits." << std::endl;
    const auto[traits, phegen_genes] = loadPheGenIGenesAndTraits(path_file_phegeni, genes);

    std::cout << "Creating gene modules." << std::endl;
    const modules gene_modules = loadPheGenIGeneModules(path_file_phegeni, genes, traits, path_file_gene_interactions);

    std::cout << "Creating protein modules." << std::endl;
    bidirectional_mapping mapping_genes_to_proteins = readMapping(path_file_mapping_proteins_to_genes, true, false);;
    const modules protein_modules = createPheGenIModules(gene_modules, genes, proteins, traits,
                                                         mapping_genes_to_proteins.second_to_first,
                                                         path_file_protein_interactions);

    std::cout << "Creating proteoform modules." << std::endl;
    bidirectional_mapping mapping_proteins_to_proteoforms = readMapping(path_file_mapping_proteins_to_proteoforms, true,
                                                                        true);
    const modules proteoform_modules = createPheGenIModules(protein_modules, proteins, proteoforms, traits,
                                                            mapping_proteins_to_proteoforms.second_to_first,
                                                            path_file_proteoform_interactions);

    writeModules(static_cast<std::string>(path_modules) + "gene_modules.tsv", gene_modules, traits, genes);
    writeModules(static_cast<std::string>(path_modules) + "protein_modules.tsv", protein_modules, traits, proteins);
    writeModules(static_cast<std::string>(path_scores) + "proteoform_modules.tsv", proteoform_modules, traits,
                 proteoforms);

    // Calculate overlap scores between all gene module pairs
    std::cout << "Calculating overlap similarity for Genes." << std::endl;
    auto gene_scores = getScores(gene_modules.group_to_members, getOverlapSimilarity);
    writeScores(gene_modules.group_to_members, gene_scores,
                path_scores + "scores_genes_overlap_similarity.tsv"));
    std::cout << "Calculating Jaccard similarity for Genes." << std::endl;
    writeScores(gene_modules.group_to_members, getScores(gene_modules.group_to_members, getJaccardSimilarity),
                path_scores + "scores_genes_jaccard_similarity.tsv"));


    // Calculate overlap scores between all protein module pairs
    std::cout << "Calculating overlap similarity for Proteins." << std::endl;
    auto protein_scores = getScores(protein_modules.group_to_members, getOverlapSimilarity);
    writeScores(protein_modules.group_to_members, protein_scores,
                path_scores + "scores_proteins_overlap_similarity.tsv"));
    std::cout << "Calculating Jaccard similarity for Proteins." << std::endl;
    writeScores(protein_modules.group_to_members, getScores(protein_modules.group_to_members, getJaccardSimilarity),
                path_scores + "scores_proteins_jaccard_similarity.tsv"));


    // Calculate overlap scores between all proteoform module pairs
    std::cout << "Calculating overlap similarity for Proteoforms." << std::endl;
    auto proteoform_scores = getScores(proteoform_modules.group_to_members, getOverlapSimilarity);
            writeScores(proteoform_modules.group_to_members, proteoform_scores,
    path_scores + "scores_proteoforms_overlap_similarity.tsv"));
    std::cout << "Calculating Jaccard similarity for Proteoforms." << std::endl;
    writeScores(proteoform_modules.group_to_members,
                getScores(proteoform_modules.group_to_members, getJaccardSimilarity),
                path_scores + "scores_proteoforms_jaccard_similarity.tsv"));

    std::cout << "Finished with the Overlap analysis";
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
