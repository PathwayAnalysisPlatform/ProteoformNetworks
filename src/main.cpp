#include <iostream>
#include <string>

#include "3_artefactual_overlap/artefactual_overlap.hpp"

using namespace std;

// Input files
const string path_file_gene_search = "resources/3_artefactual_overlap/gene_search.tsv";
const string path_file_protein_search = "resources/3_artefactual_overlap/protein_search.tsv";
const string path_file_proteoform_search = "resources/3_artefactual_overlap/proteoform_search.tsv";

// Output files
const string path_file_report_degree_reduction_analysis = "reports/1_degree_reduction_analysis.txt";
const string path_file_report_percolation_analysis = "reports/2_percolation_analysis.txt";
const string path_file_report_artefactual_overlap_analysis = "reports/3_artefactual_overlap_analysis.txt";
const string path_file_gene_art_pairs = "reports/pathway_pairs_with_gene_art_overlap.txt";
const string path_file_protein_art_pairs = "reports/pathway_pairs_with_protein_art_overlap.txt";
const string path_file_report_key_ptm_analysis = "reports/4_key_ptm_analysis.txt";
const string path_file_pathways_in_key_ptm_overlaps = "reports/4_key_ptms_overlaps/pathway_pairs.txt";
const string path_file_proteins_in_key_ptm_overlaps = "reports/4_key_ptms_overlaps/proteins.txt";
const string path_file_proteoforms_in_key_ptm_overlaps = "reports/4_key_ptms_overlaps/proteoforms.txt";
const string path_file_modifications_in_key_ptm_overlaps = "reports/4_key_ptms_overlaps/modifications.txt";

void doDegreeReductionAnalysis() {
    // TODO
}

void doPercolationAnalysis() {
    // TODO
}

void doKeyPtmAnalysis() {
    // set<pair<string, string>> pairs = findPairsWithKeyPTMExamples(
    //     0.0,
    //     path_file_proteoform_search,
    //     path_file_pathways_in_key_ptm_overlaps,
    //     path_file_proteins_in_key_ptm_overlaps,
    //     path_file_proteoforms_in_key_ptm_overlaps,
    //     path_file_modifications_in_key_ptm_overlaps);
}

int main()
try {
    doDegreeReductionAnalysis();  // 1)
    doPercolationAnalysis();      // 2)
    doArtefactualOverlapAnalysis(path_file_gene_search,
                                 path_file_protein_search,
                                 path_file_proteoform_search,
                                 path_file_report_artefactual_overlap_analysis,
                                 path_file_gene_art_pairs,
                                 path_file_protein_art_pairs);  // 3)
    doKeyPtmAnalysis();                                         // 4)
} catch (const std::exception& ex) {
   cout << ex.what() << "\n";
}
