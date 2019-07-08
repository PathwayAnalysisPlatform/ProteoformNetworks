#ifndef REACTOME_HPP_
#define REACTOME_HPP_

#include <unordered_map>
#include <bitset>
#include <cstdio>

#include "types.hpp"
#include "bimap.hpp"
#include "entity.hpp"

// Reactome v69
const size_t REACTOME_GENES = 24109;
const size_t REACTOME_PROTEINS = 10833;
const size_t REACTOME_PROTEOFORMS = 13991;

bimap loadReactomeEntities(std::string_view path_search_file);
vs loadReactomeEntitiesIndexToEntitites(std::string_view path_search_file);

struct load_mapping_proteins_proteoforms_result {
	ummss protein_to_proteoforms;
	ummss proteoform_to_proteins;
};

load_mapping_proteins_proteoforms_result loadMappingProteinsProteoforms(std::string_view path_file_mapping);

using reactome_gene_sets = um<std::string, std::bitset<REACTOME_GENES>>;
using reactome_protein_sets = um<std::string, std::bitset<REACTOME_PROTEINS>>;
using reactome_proteoform_sets = um<std::string, std::bitset<REACTOME_PROTEOFORMS>>;

reactome_gene_sets loadGeneSets(std::string_view file_path, const bimap& genes, bool pathways);
reactome_protein_sets loadProteinSets(std::string_view file_path, const bimap& proteins, bool pathways);
reactome_proteoform_sets loadProteoformSets(std::string_view file_path, const bimap& proteoforms, bool pathways);

umss loadPathwayNames(std::string_view path_search_file);

bimap deductProteinsFromGenes(
	std::string_view path_file_mapping_proteins_genes,
	const ummss& genes_to_proteins,
	const bimap& genes);
bimap deductProteoformsFromProteins(
	const ummss& proteins_to_proteoforms, 
	const bimap& proteins);

struct reactome_networks {
	ummss adjacency_list_proteins;
	ummss adjacency_list_proteoforms;
};

reactome_networks loadReactomeNetworks(
	std::string_view path_file_gene_search,
	std::string_view path_file_protein_search,
	std::string_view path_file_proteoform_search);

ummss loadProteinsAdjacencyList(std::string_view search_file_path);

ummss loadProteoformsAdjacencyList(std::string_view search_file_path);

#endif // !REACTOME_HPP_

