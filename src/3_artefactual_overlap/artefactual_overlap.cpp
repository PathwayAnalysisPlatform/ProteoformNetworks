#include "artefactual_overlap.hpp"

using namespace std;

namespace artefactual_overlap {

/**
 * Find artefactual overlaps: pairs of pathways that share nodes only in the
 * gene or protein level, but not at the proteoform level
 */
set<pair<string, string>> findPathwayPairs(const map<pair<string, string>, bitset<NUM_GENES>>& overlapping_gene_set_pairs,
                                           const map<pair<string, string>, bitset<NUM_PROTEINS>>& overlapping_protein_set_pairs,
                                           const map<pair<string, string>, bitset<NUM_PROTEOFORMS>>& overlapping_proteoform_set_pairs) {
   set<pair<string, string>> result;

   auto t0 = clock();
   cerr << "Comparing gene and proteoform pairs..." << endl;
   for (const auto& gene_set_pair : overlapping_gene_set_pairs) {
      auto it = overlapping_proteoform_set_pairs.find(gene_set_pair.first);
      if (it == overlapping_proteoform_set_pairs.end()) {
         result.insert(gene_set_pair.first);
      }
   }
   auto t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";

   t0 = clock();
   cerr << "Comparing protein and proteoform pairs..." << endl;
   for (const auto& protein_set_pair : overlapping_protein_set_pairs) {
      auto it = overlapping_proteoform_set_pairs.find(protein_set_pair.first);
      if (it == overlapping_proteoform_set_pairs.end()) {
         result.insert(protein_set_pair.first);
      }
   }
   t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";

   return result;
}

bitset<NUM_PROTEOFORMS> getProteoformsWithAccessions(const set<string>& accessions,
                                                     const bitset<NUM_PROTEOFORMS>& proteoform_set,
                                                     const vector<string>& index_to_proteoforms) {
   bitset<NUM_PROTEOFORMS> result;
   for (int I = 0; I < proteoform_set.size(); I++) {
      if (proteoform_set.test(I)) {
         if (accessions.find(getAccession(index_to_proteoforms.at(I))) != accessions.end()) {
            result.set(I);
         }
      }
   }
   return result;
}

void writeReportRecords(ofstream& output,
                        const set<pair<string, string>> examples,
                        const map<string, string>& sets_to_names,
                        const map<string, bitset<NUM_GENES>>& sets_to_genes,
                        const map<string, bitset<NUM_PROTEINS>>& sets_to_proteins,
                        const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_proteoforms,
                        const vector<string> index_to_genes,
                        const vector<string> index_to_proteins,
                        const vector<string> index_to_proteoforms) {
   auto t0 = clock();

   for (const auto& example : examples) {
      bitset<NUM_GENES> overlap_genes = sets_to_genes.at(example.first) & sets_to_genes.at(example.second);
      bitset<NUM_PROTEINS> overlap_proteins = sets_to_proteins.at(example.first) & sets_to_proteins.at(example.second);
      bitset<NUM_PROTEOFORMS> overlap_proteoforms = sets_to_proteoforms.at(example.first) & sets_to_proteoforms.at(example.second);
      auto decomposed_overlap_proteoforms_1 = getProteoformsWithAccessions(getProteinStrings(overlap_proteins, index_to_proteins),
                                                                           sets_to_proteoforms.at(example.first), index_to_proteoforms);
      auto decomposed_overlap_proteoforms_2 = getProteoformsWithAccessions(getProteinStrings(overlap_proteins, index_to_proteins),
                                                                           sets_to_proteoforms.at(example.second), index_to_proteoforms);

      output << example.first << "\t" << example.second << "\t" << sets_to_names.at(example.first) << "\t" << sets_to_names.at(example.second) << "\t";
      output << sets_to_genes.at(example.first).count() << "\t" << sets_to_proteins.at(example.first).count() << "\t";
      output << sets_to_proteoforms.at(example.first).count() << "\t";
      output << sets_to_genes.at(example.second).count() << "\t" << sets_to_proteins.at(example.second).count() << "\t";
      output << sets_to_proteoforms.at(example.second).count() << "\t";
      output << overlap_genes.count() << "\t" << overlap_proteins.count() << "\t" << overlap_proteoforms.count() << "\t";
      printGeneMembers(output, overlap_genes, index_to_genes);
      output << "\t";
      printProteinMembers(output, overlap_proteins, index_to_proteins);
      output << "\t";
      printProteoformMembers(output, decomposed_overlap_proteoforms_1, index_to_proteoforms);
      output << "\t";
      printProteoformMembers(output, decomposed_overlap_proteoforms_2, index_to_proteoforms);
      output << "\n";
   }
   auto t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
}

void writePathwayReport(ofstream& output,
                        const set<pair<string, string>> examples,
                        const map<string, string>& sets_to_names,
                        const map<string, bitset<NUM_GENES>>& sets_to_genes,
                        const map<string, bitset<NUM_PROTEINS>>& sets_to_proteins,
                        const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_proteoforms,
                        const vector<string> index_to_genes,
                        const vector<string> index_to_proteins,
                        const vector<string> index_to_proteoforms) {
   cerr << "Reporting pathway pairs with artefactual overlap...\n";
   output << "PATHWAY_1\tPATHWAY_2\tPATHWAY_1_NAME\tPATHWAY_2_NAME\t";
   output << "PATHWAY_1_GENE_SIZE\tPATHWAY_1_PROTEIN_SIZE\tPATHWAY_1_PROTEOFORM_SIZE\t";
   output << "PATHWAY_2_GENE_SIZE\tPATHWAY_2_PROTEIN_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\t";
   output << "GENE_OVERLAP\tPROTEIN_OVERLAP\tPROTEOFORM_OVERLAP\t";
   output << "OVERLAP_GENES\tOVERLAP_PROTEINS\tDECOMPOSED_OVERLAP_PROTEOFORMS_1\tDECOMPOSED_OVERLAP_PROTEOFORMS_2\n";
   writeReportRecords(output, examples, sets_to_names,
                      sets_to_genes, sets_to_proteins, sets_to_proteoforms,
                      index_to_genes, index_to_proteins, index_to_proteoforms);
}

void writePhenotypeReport(ofstream& output,
                          const set<pair<string, string>> examples,
                          const map<string, string>& sets_to_names,
                          const map<string, bitset<NUM_GENES>>& sets_to_genes,
                          const map<string, bitset<NUM_PROTEINS>>& sets_to_proteins,
                          const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_proteoforms,
                          const vector<string> index_to_genes,
                          const vector<string> index_to_proteins,
                          const vector<string> index_to_proteoforms) {
   output << "PHENOTYPE_1\tPHENOTYPE_2\tPHENOTYPE_1_NAME\tPHENOTYPE_2_NAME\t";
   output << "PHENOTYPE_1_GENE_SIZE\tPHENOTYPE_1_PROTEIN_SIZE\tPHENOTYPE_1_PROTEOFORM_SIZE\t";
   output << "PHENOTYPE_2_GENE_SIZE\tPHENOTYPE_2_PROTEIN_SIZE\tPHENOTYPE_2_PROTEOFORM_SIZE\t";
   output << "GENE_OVERLAP\tPROTEIN_OVERLAP\tPROTEOFORM_OVERLAP\t";
   output << "OVERLAP_GENES\tOVERLAP_PROTEINS\tDECOMPOSED_OVERLAP_PROTEOFORMS_1\tDECOMPOSED_OVERLAP_PROTEOFORMS_2\n";
   writeReportRecords(output, examples, sets_to_names,
                      sets_to_genes, sets_to_proteins, sets_to_proteoforms,
                      index_to_genes, index_to_proteins, index_to_proteoforms);
}

void reportPathwayPairs(const string& path_file_gene_search,
                        const string& path_file_protein_search,
                        const string& path_file_proteoform_search,
                        const string& report_file_path) {
   const auto [index_to_genes, genes_to_index] = loadEntities(path_file_gene_search);
   const auto [index_to_proteins, proteins_to_index] = loadEntities(path_file_protein_search);
   const auto [index_to_proteoforms, proteoforms_to_index] = loadEntities(path_file_proteoform_search);
   const auto pathways_to_names = loadPathwayNames(path_file_proteoform_search);

   map<string, bitset<NUM_GENES>> pathways_to_genes = loadPathwaysGeneMembers(path_file_gene_search, genes_to_index);
   map<string, bitset<NUM_PROTEINS>> pathways_to_proteins = loadPathwaysProteinMembers(path_file_protein_search, proteins_to_index);
   map<string, bitset<NUM_PROTEOFORMS>> pathways_to_proteoforms = loadPathwaysProteoformMembers(path_file_proteoform_search, proteoforms_to_index);

   cout << "Calculating gene sets overlap..." << endl;
   const auto overlapping_gene_set_pairs = findOverlappingGeneSets(pathways_to_genes, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_PATHWAY_SIZE, MAX_PATHWAY_SIZE);
   cout << "Calculating protein sets overlap..." << endl;
   const auto overlapping_protein_set_pairs = findOverlappingProteinSets(pathways_to_proteins, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_PATHWAY_SIZE, MAX_PATHWAY_SIZE);
   cout << "Calculating proteoform sets..." << endl;
   const auto overlapping_proteoform_set_pairs = findOverlappingProteoformSets(pathways_to_proteoforms);

   cout << "Finding examples of artifactual overlap...\n";
   const auto examples = findPathwayPairs(overlapping_gene_set_pairs,
                                          overlapping_protein_set_pairs,
                                          overlapping_proteoform_set_pairs);

   ofstream report(report_file_path);
   writePathwayReport(report, examples, pathways_to_names,
                      pathways_to_genes, pathways_to_proteins, pathways_to_proteoforms,
                      index_to_genes, index_to_proteins, index_to_proteoforms);
}

void reportPhenotypePairs() {
   // TODO: Report pairs of disease modules
}

// Find pairs of modules/pathways that overlap on gene or protein network, but not in proteoform network
void doAnalysis(const string& path_file_gene_search,
                const string& path_file_protein_search,
                const string& path_file_proteoform_search,
                const string& report_file_path) {
   cout << "Searching for artefactual overlap examples..." << endl;

   reportPathwayPairs(path_file_gene_search,
                      path_file_protein_search,
                      path_file_proteoform_search,
                      report_file_path);

   reportPhenotypePairs();
}

}  // namespace artefactual_overlap
