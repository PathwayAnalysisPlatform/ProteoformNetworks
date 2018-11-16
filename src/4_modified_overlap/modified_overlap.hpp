#ifndef MODIFIED_OVERLAP_H_
#define MODIFIED_OVERLAP_H_

#include "../overlap.hpp"

#include <regex>
#include <charconv>

const float MIN_MODIFIED_ALL_MEMBERS_RATIO = 0.1;
const float MIN_MODIFIED_OVERLAP_MEMBERS_RATIO = 0.9;

const std::regex RGX_MODIFICATION{"[;,]\\d{5}"};

namespace modified_overlap {
void doAnalysis(const std::string& path_file_proteoform_search,
                const std::string& path_file_PheGenI_full,
                const std::string& path_file_mapping_proteins_to_genes,
                const std::string& path_file_report_pathway,
                const std::string& path_file_modified_overlap_pathway_proteins,
                const std::string& path_file_modified_overlap_pathway_proteoforms,
                const std::string& path_file_modified_overlap_pathway_modifications,
                const std::string& path_file_report_trait,
                const std::string& path_file_modified_overlap_trait_proteins,
                const std::string& path_file_modified_overlap_trait_proteoforms,
                const std::string& path_file_modified_overlap_trait_modifications);
}

#endif /* MODIFIED_OVERLAP_H_ */