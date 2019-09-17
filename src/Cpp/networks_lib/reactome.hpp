#ifndef REACTOME_HPP_
#define REACTOME_HPP_

#include <bitset>
#include <cstring>

#include "types.hpp"
#include "bimap_str_int.hpp"

entity_mapping loadMappingProteinsProteoforms(std::string_view path_file_mapping);

umsb loadGeneSets(std::string_view file_path, const bimap_str_int &phegeni_genes, bool pathways);

umsb loadProteinSets(std::string_view file_path, const bimap_str_int &proteins, bool pathways);

umsb loadProteoformSets(std::string_view file_path, const bimap_str_int &proteoforms, bool pathways);

umss loadPathwayNames(std::string_view path_search_file);

#endif // !REACTOME_HPP_

