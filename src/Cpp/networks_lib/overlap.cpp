#include "overlap.hpp"

using namespace std;

void doOverlapAnalysis(
        std::string_view path_file_PheGenI,
        std::string_view path_file_reactome_genes,
        std::string_view path_file_mapping_proteins_to_genes,
        std::string_view path_file_protein_search,
        std::string_view path_file_proteoform_search) {

    // Read Reactome genes. Take them as all acceptable gene names.
    const bimap_str_int reactome_genes = createBimap(path_file_reactome_genes); // Gene names --> str_to_int

    // Read Phegeni trait modules with genes as members. Only gene members also in the acceptable gene list.
    const auto modules = loadPheGenIModules(path_file_PheGenI, reactome_genes);

    // Calculate overlap scores between all trait gene set pairs
    const auto scores_overlap_similarity = getScores(modules.traits_to_genes, getOverlapSimilarity);
    const auto scores_jaccard_similarity = getScores(modules.traits_to_genes, getJaccardSimilarity);



    // Create proteoform sets for each trait
//	const auto [adjacency_list_proteins, adjacency_list_proteoforms] = loadReactomeNetworks(path_file_protein_search, path_file_proteoform_search);
//
//	const auto [genes_to_proteins, proteins_to_genes] = loadMappingGenesProteins(path_file_mapping_proteins_to_genes);
//	const auto phegeni_proteins = deductProteinsFromGenes(genes_to_proteins, phegeni_genes);
//	const um<string, bitset<PHEGENI_PROTEINS>> traits_to_proteins = convertGeneSets(modules.traits_to_genes, phegeni_genes, genes_to_proteins, phegeni_proteins, adjacency_list_proteins);
//
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
