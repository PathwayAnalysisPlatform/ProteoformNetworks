#include "overlap_analysis.hpp"

using namespace std;

void doOverlapAnalysis(
        std::string_view path_file_phegeni,
        std::string_view path_file_reactome_genes,
        std::string_view path_file_reactome_proteins,
        std::string_view path_file_reactome_proteoforms,
        std::string_view path_file_mapping_proteins_to_genes,
        std::string_view path_file_mapping_proteins_to_proteoforms,
        std::string_view path_file_protein_interactions,
        std::string_view path_file_proteoform_interactions,
        std::string_view path_scores) {

    // Read Reactome genes, proteins and proteoforms. Take them as set of acceptable identifiers.
    const bimap_str_int genes = createBimap(path_file_reactome_genes);
    const bimap_str_int proteins = createBimap(path_file_reactome_proteins);
    const bimap_str_int proteoforms = createBimap(path_file_reactome_proteoforms);

    // Read traits and genes in Phegeni file
    const auto[phegeni_genes, traits] = loadPheGenIGenesAndTraits(path_file_phegeni, genes);

    const modules gene_modules = loadPheGenIGeneModules(path_file_phegeni, genes, traits);

    bidirectional_mapping mapping_genes_to_proteins = readMapping(path_file_mapping_proteins_to_genes, true, false);;
    const modules protein_modules = createPheGenIModules(gene_modules, genes, proteins, traits,
                                                         mapping_genes_to_proteins.second_to_first,
                                                         path_file_protein_interactions);

    bidirectional_mapping mapping_proteins_to_proteoforms = readMapping(path_file_mapping_proteins_to_proteoforms, true,
                                                                        true);
    const modules proteoform_modules = createPheGenIModules(protein_modules, proteins, proteoforms, traits,
                                                            mapping_proteins_to_proteoforms.first_to_second,
                                                            path_file_proteoform_interactions);

    // Calculate overlap scores between all gene module pairs
    std::string file_name = "scores_genes_overlap_similarity.tsv";
    file_name = path_scores.data() + file_name;
    writeScores(gene_modules.group_to_members, getScores(gene_modules.group_to_members, getOverlapSimilarity),
                file_name);
    file_name = "scores_genes_jaccard_similarity.tsv";
    file_name = path_scores.data() + file_name;
    writeScores(gene_modules.group_to_members, getScores(gene_modules.group_to_members, getJaccardSimilarity),
                file_name);

    // Calculate overlap scores between all protein module pairs
    file_name = "scores_proteins_overlap_similarity.tsv";
    file_name = path_scores.data() + file_name;
    writeScores(gene_modules.group_to_members, getScores(gene_modules.group_to_members, getOverlapSimilarity),
                file_name);
    file_name = "scores_proteins_jaccard_similarity.tsv";
    file_name = path_scores.data() + file_name;
    writeScores(gene_modules.group_to_members, getScores(gene_modules.group_to_members, getJaccardSimilarity),
                file_name);

    // Calculate overlap scores between all proteoform module pairs
    file_name = "scores_proteoforms_overlap_similarity.tsv";
    file_name = path_scores.data() + file_name;
    writeScores(gene_modules.group_to_members, getScores(gene_modules.group_to_members, getOverlapSimilarity),
                file_name);
    file_name = "scores_proteoforms_jaccard_similarity.tsv";
    file_name = path_scores.data() + file_name;
    writeScores(gene_modules.group_to_members, getScores(gene_modules.group_to_members, getJaccardSimilarity),
                file_name);

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
