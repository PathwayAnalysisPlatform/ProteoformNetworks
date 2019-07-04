#ifndef MODIFIED_OVERLAP_H_
#define MODIFIED_OVERLAP_H_

#include <charconv>
#include <regex>

#include "overlap.hpp"
#include "proteoform.hpp"

const double MIN_MODIFIED_ALL_MEMBERS_RATIO = 0.1;
const double MIN_MODIFIED_OVERLAP_MEMBERS_RATIO = 0.9;

namespace modified_overlap {
void doAnalysis(std::string_view path_file_gene_search,
                std::string_view path_file_protein_search,
                std::string_view path_file_proteoform_search,
                std::string_view path_file_PheGenI,
                std::string_view path_file_mapping_proteins_to_genes,
                std::string_view path_file_report_pathway,
                std::string_view path_file_modified_overlap_pathway_proteins,
                std::string_view path_file_modified_overlap_pathway_proteoforms,
                std::string_view path_file_modified_overlap_pathway_modifications,
                std::string_view path_file_report_trait,
                std::string_view path_file_modified_overlap_trait_proteins,
                std::string_view path_file_modified_overlap_trait_proteoforms,
                std::string_view path_file_modified_overlap_trait_modifications);
}

#endif /* MODIFIED_OVERLAP_H_ */