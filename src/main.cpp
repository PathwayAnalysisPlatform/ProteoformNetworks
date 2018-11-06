#include "1_degree_reduction/degree_reduction.hpp"
#include "2_percolation/percolation.hpp"
#include "3_artefactual_overlap/artefactual_overlap.hpp"
#include "4_modified_overlap/modified_overlap.hpp"

// Input files
const std::string path_file_gene_search = "resources/3_artefactual_overlap/gene_search.tsv";
const std::string path_file_protein_search = "resources/3_artefactual_overlap/protein_search.tsv";
const std::string path_file_proteoform_search = "resources/3_artefactual_overlap/proteoform_search.tsv";

// Output files
const std::string path_file_report_degree_reduction_analysis = "reports/1_degree_reduction_analysis.txt";
const std::string path_file_report_percolation_analysis = "reports/2_percolation_analysis.txt";
const std::string path_file_report_artefactual_overlap_analysis = "reports/3_artefactual_overlap_analysis.txt";
const std::string path_file_report_modified_overlap_analysis = "reports/4_modified_overlap_analysis.txt";
const std::string path_file_proteins_in_modified_overlap = "reports/4_modified_overlap/proteins.txt";
const std::string path_file_proteoforms_in_modified_overlap = "reports/4_modified_overlap/proteoforms.txt";
const std::string path_file_modifications_in_modified_overlap = "reports/4_modified_overlap/modifications.txt";

int main() try {
   doDegreeReductionAnalysis();
   doPercolationAnalysis();
   doArtefactualOverlapAnalysis(path_file_gene_search,
                                path_file_protein_search,
                                path_file_proteoform_search,
                                path_file_report_artefactual_overlap_analysis);
   doModifiedOverlapAnalysis(path_file_proteoform_search,
                             path_file_report_modified_overlap_analysis,
                             path_file_proteins_in_modified_overlap,
                             path_file_proteoforms_in_modified_overlap,
                             path_file_modifications_in_modified_overlap);
} catch (const std::exception& ex) {
   std::cout << ex.what() << "\n";
}