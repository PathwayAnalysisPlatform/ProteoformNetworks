#include "artefactual_overlap.hpp"
#include "degree_reduction.hpp"
#include "lib/pathway/dataset.hpp"
#include "modified_overlap.hpp"
#include "percolation.hpp"

// ************* Input files

// 1) Degree reduction
const std::string path_file_num_pathways_per_protein = "resources/reactome/num_pathways_per_protein.csv";
const std::string path_file_num_pathways_per_proteoform = "resources/reactome/num_pathways_per_proteoform.csv";
const std::string path_file_gene_mapping = "resources/reactome/gene_mapping.csv";              // Includes the mapping from genes to proteins and genes to reactions and pathways.
const std::string path_file_protein_mapping = "resources/reactome/protein_mapping.csv";        // Includes the mapping from proteins to reactions and pathways.
const std::string path_file_proteoform_mapping = "resources/reactome/proteoform_mapping.csv";  // Includes the mapping from proteoforms to reactions and pathways.

// 2) Percolation

// 3) Artefactual overlap
const std::string path_file_gene_search = "resources/3_artefactual_overlap/gene_search.tsv";
const std::string path_file_protein_search = "resources/3_artefactual_overlap/protein_search.tsv";
const std::string path_file_proteoform_search = "resources/3_artefactual_overlap/proteoform_search.tsv";
const std::string path_file_PheGenI_full = "resources/PheGenI/PheGenI_Association_full.tab";
const std::string path_file_mapping_proteins_to_genes = "resources/UniProt/proteins_to_genes.tab";

// 4) Modified overlap

// ************* Output files
const std::string path_file_report_degree_reduction_analysis = "reports/degree_reduction_analysis.txt";
const std::string path_file_report_percolation_analysis = "reports/percolation_analysis.txt";

const std::string path_file_report_artefactual_overlap_pathway = "reports/artefactual_overlap_pathway.txt";
const std::string path_file_report_artefactual_overlap_trait = "reports/artefactual_overlap_trait.txt";

const std::string path_file_report_modified_overlap_pathway = "reports/modified_overlap_pathway.txt";
const std::string path_file_modified_overlap_pathway_proteins = "reports/modified_overlap_pathway_proteins.txt";
const std::string path_file_modified_overlap_pathway_proteoforms = "reports/modified_overlap_pathway_proteoforms.txt";
const std::string path_file_modified_overlap_pathway_modifications = "reports/modified_overlap_pathway_modifications.txt";

const std::string path_file_report_modified_overlap_trait = "reports/modified_overlap_trait.txt";
const std::string path_file_modified_overlap_trait_proteins = "reports/modified_overlap_trait_proteins.txt";
const std::string path_file_modified_overlap_trait_proteoforms = "reports/modified_overlap_trait_proteoforms.txt";
const std::string path_file_modified_overlap_trait_modifications = "reports/modified_overlap_trait_modifications.txt";

int main() try {
   pathway::dataset pathwayDataSet(path_file_gene_mapping, path_file_protein_mapping, path_file_proteoform_mapping);

   degree_reduction::doAnalysis(pathwayDataSet, path_file_report_degree_reduction_analysis);
   std::cout << "\n\n-----------------------------------------------------\n\n";
   // percolation::doAnalysis(path_file_gene_search,
   //                         path_file_protein_search,
   //                         path_file_proteoform_search);
   std::cout << "\n\n-----------------------------------------------------\n\n";
   artefactual_overlap::doAnalysis(path_file_gene_search,
                                   path_file_protein_search,
                                   path_file_proteoform_search,
                                   path_file_PheGenI_full,
                                   path_file_mapping_proteins_to_genes,
                                   path_file_report_artefactual_overlap_pathway,
                                   path_file_report_artefactual_overlap_trait);
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
   system("pause");
} catch (const std::exception& ex) {
   std::cout << ex.what() << "\n";
}