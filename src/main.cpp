#include "1_degree_reduction/degree_reduction.hpp"
#include "2_percolation/percolation.hpp"
#include "3_rule_out_gene_centric_overlap/rule_out_gene_centric_overlap.hpp"
#include "4_modified_overlap/modified_overlap.hpp"

// Input files
const std::string path_file_gene_search = "resources/3_rule_out_gene_centric_overlap/gene_search.tsv";
const std::string path_file_protein_search = "resources/3_rule_out_gene_centric_overlap/protein_search.tsv";
const std::string path_file_proteoform_search = "resources/3_rule_out_gene_centric_overlap/proteoform_search.tsv";

const std::string path_file_PheGenI_full = "resources/PheGenI/PheGenI_Association_full.tab";
const std::string path_file_mapping_proteins_to_genes = "resources/UniProt/proteins_to_genes.tab";

// Output files
const std::string path_file_report_degree_reduction_analysis = "reports/degree_reduction_analysis.txt";
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

int main() try {
   doDegreeReductionAnalysis();

   doPercolationAnalysis();

   rule_out_gene_centric_overlap::doAnalysis(path_file_gene_search,
                                   path_file_protein_search,
                                   path_file_proteoform_search,
                                   path_file_PheGenI_full,
                                   path_file_mapping_proteins_to_genes,
                                   path_file_report_rule_out_gene_centric_overlap_pathway,
                                   path_file_report_rule_out_gene_centric_overlap_trait);
   std::cout << "\n\n-----------------------------------------------------\n\n";
   modified_overlap::doAnalysis(path_file_gene_search,
                                path_file_protein_search,
                                path_file_proteoform_search,
                                path_file_PheGenI_full,
                                path_file_mapping_proteins_to_genes,
                                path_file_report_modified_overlap_pathway,
                                path_file_modified_overlap_pathway_proteins,
                                path_file_modified_overlap_pathway_proteoforms,
                                path_file_modified_overlap_pathway_modifications,
                                path_file_report_modified_overlap_trait,
                                path_file_modified_overlap_trait_proteins,
                                path_file_modified_overlap_trait_proteoforms,
                                path_file_modified_overlap_trait_modifications);
   // system("pause");
} catch (const std::exception& ex) {
   std::cout << ex.what() << "\n";
}