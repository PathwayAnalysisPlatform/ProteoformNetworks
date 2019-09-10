#ifndef UNIPROT_HPP_
#define UNIPROT_HPP_

#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string_view>

#include "types.hpp"

struct load_mapping_genes_proteins_result {
	ummss gene_to_proteins;
	ummss protein_to_genes;
};

load_mapping_genes_proteins_result loadMappingGenesProteins(std::string_view path_file_mapping);

#endif // !UNIPROT_HPP_

