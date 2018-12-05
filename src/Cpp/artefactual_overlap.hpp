#ifndef ARTEFACTUAL_OVERLAP_H_
#define ARTEFACTUAL_OVERLAP_H_

#include <string_view>

#include "overlap.hpp"
#include "lib/proteoform/proteoform.h"

namespace artefactual_overlap {
void doAnalysis(std::string_view path_file_gene_search,
                std::string_view path_file_protein_search,
                std::string_view path_file_proteoform_search,
                std::string_view path_file_PheGenI_full,
                std::string_view path_file_mapping_proteins_to_genes,
                std::string_view path_file_report_pathway,
                std::string_view path_file_report_trait);
}

#endif /* ARTEFACTUAL_OVERLAP_H_ */