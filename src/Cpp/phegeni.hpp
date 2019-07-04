#ifndef PHEGENI_HPP
#define PHEGENI_HPP

#include <fstream>
#include <string>

#include "bimap.hpp"

const long double GENOME_WIDE_SIGNIFICANCE = 5e-8;

struct load_phegeni_entities_result {
	bimap genes;
	bimap traits;
};

struct load_phegeni_entities_index_to_entitites_result {
	vs index_to_genes;
	vs index_to_traits;
};

load_phegeni_entities_result loadPheGenIEntities(std::string_view path_file_PheGenI,
	const double& max_p_value,
	const umsi& reactome_genes_to_index);

load_phegeni_entities_index_to_entitites_result loadPheGenIEntitiesIndexToEntitites(const std::string& path_file_PheGenI,
	const double& max_p_value,
	const umsi& reactome_genes_to_index);

#endif // !PHEGENI_HPP

