#ifndef GENE_LEVEL_ONLY_OVERLAP_H
#define GENE_LEVEL_ONLY_OVERLAP_H

#include <string_view>

#include "overlap.hpp"
#include "proteoform.hpp"
#include "reactome.hpp"
#include "phegeni.hpp"

namespace gene_level_only_overlap {
void doAnalysis(std::string_view path_file_gene_search,
                std::string_view path_file_protein_search,
                std::string_view path_file_proteoform_search,
                std::string_view path_file_PheGenI,
                std::string_view path_file_mapping_proteins_to_genes,
                std::string_view path_file_report_pathway,
                std::string_view path_file_report_trait);
}

#endif /* GENE_LEVEL_ONLY_OVERLAP_H */