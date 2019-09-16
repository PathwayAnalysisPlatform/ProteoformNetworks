#ifndef REACTOME_HPP_
#define REACTOME_HPP_

#include <unordered_map>
#include <bitset>
#include <cstring>

#include "types.hpp"
#include "bimap_str_int.hpp"

entity_mapping loadMappingProteinsProteoforms(std::string_view path_file_mapping);

umsb loadGeneSets(std::string_view file_path, const bimap_str_int& phegeni_genes, bool pathways);
umsb loadProteinSets(std::string_view file_path, const bimap_str_int& proteins, bool pathways);
umsb loadProteoformSets(std::string_view file_path, const bimap_str_int& proteoforms, bool pathways);

umss loadPathwayNames(std::string_view path_search_file);

bimap_str_int deductProteinsFromGenes(
	const ummss& genes_to_proteins,
	const bimap_str_int& phegeni_genes);

bimap_str_int deductProteoformsFromProteins(
	const ummss& proteins_to_proteoforms, 
	const bimap_str_int& proteins);

#endif // !REACTOME_HPP_

