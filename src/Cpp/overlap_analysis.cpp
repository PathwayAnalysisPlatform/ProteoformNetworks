#include "overlap_analysis.hpp"

using namespace std;

void doOverlapAnalysis(
        std::string_view path_file_phegeni,
        std::string_view path_file_reactome_genes,
        std::string_view path_file_reactome_proteins,
        std::string_view path_file_mapping_proteins_to_genes,
        std::string_view path_file_protein_edges,
        std::string_view path_file_proteoform_edges,
        std::string_view path_scores) {

    // Read Reactome genes. Take them as all acceptable gene names.
    const bimap_str_int genes = createBimap(path_file_reactome_genes); // Gene names --> str_to_int

    // Read traits and genes in Phegeni file
    const auto[phegeni_genes, traits] = loadPheGenIGenesAndTraits(path_file_phegeni, genes);

    // Read Phegeni trait modules with genes as members. Only gene members also in the acceptable gene list.
    const auto gene_modules = loadPheGenIGeneModules(path_file_phegeni, genes, traits);

    // Calculate overlap scores between all trait gene set pairs
    std::string file_name = "scores_overlap_similarity.tsv";
    file_name = path_scores.data() + file_name;
    writeScores(gene_modules.traits_to_entities, getScores(gene_modules.traits_to_entities, getOverlapSimilarity),
                file_name);
    file_name = "scores_jaccard_similarity.tsv";
    file_name = path_scores.data() + file_name;
    writeScores(gene_modules.traits_to_entities, getScores(gene_modules.traits_to_entities, getJaccardSimilarity),
                file_name);


    // Create proteoform sets for each trait
    const auto proteins = createBimap(path_file_reactome_proteins);
    const auto protein_modules = createPheGenIProteinModules(gene_modules,
                                                             genes,
                                                             proteins,
                                                             traits,
                                                             path_file_mapping_proteins_to_genes,
                                                             path_file_protein_edges);

//	const auto [proteins_to_proteoforms, proteoforms_to_proteins] = loadMappingProteinsProteoforms(path_file_proteoform_search);
//	const auto phegeni_proteoforms = deductProteoformsFromProteins(proteins_to_proteoforms, phegeni_proteins);
//	const um<string, bitset<PHEGENI_PROTEOFORMS>> traits_to_proteoforms = convertProteinSets(traits_to_proteins, phegeni_proteins, proteins_to_proteoforms, phegeni_proteoforms, adjacency_list_proteoforms);
//

    // Calculate overlap score between all trait proteoform set pairs

    // Write scores to file
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
