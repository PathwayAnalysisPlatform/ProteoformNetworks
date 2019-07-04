#ifndef PERCOLATION_H_
#define PERCOLATION_H_

#include "overlap.hpp"

const std::string path_file_report_percolation_analysis = "reports/percolation_analysis.txt";

const std::string path_file_report_rule_out_gene_centric_overlap_pathway = "reports/rule_out_gene_centric_overlap_pathway.txt";
const std::string path_file_report_rule_out_gene_centric_overlap_trait = "reports/rule_out_gene_centric_overlap_trait.txt";

const std::string path_file_report_modified_overlap_pathway = "reports/modified_overlap_pathway.txt";
const std::string path_file_modified_overlap_pathway_proteins = "reports/modified_overlap_pathway_proteins.txt";
const std::string path_file_modified_overlap_pathway_proteoforms = "reports/modified_overlap_pathway_proteoforms.txt";
const std::string path_file_modified_overlap_pathway_modifications = "reports/modified_overlap_pathway_modifications.txt";

const std::string path_file_report_modified_overlap_trait = "reports/modified_overlap_trait.txt";
const std::string path_file_modified_overlap_trait_proteins = "reports/modified_overlap_trait_proteins.txt";
const std::string path_file_modified_overlap_trait_proteoforms = "reports/modified_overlap_trait_proteoforms.txt";
const std::string path_file_modified_overlap_trait_modifications = "reports/modified_overlap_trait_modifications.txt";

namespace percolation {

void doAnalysis(std::string_view  path_file_gene_search,
                std::string_view  path_file_protein_search,
                std::string_view  path_file_proteoform_search);
}
#endif /* PERCOLATION_H_ */