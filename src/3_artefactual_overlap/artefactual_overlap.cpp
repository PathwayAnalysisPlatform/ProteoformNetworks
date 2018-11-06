#include "artefactual_overlap.hpp"

using namespace std;

/**
 * Find artefactual overlaps: pairs of pathways that share nodes only in the
 * gene or protein level, but not at the proteoform level
 */
set<pair<string, string>> findPathwayPairsWithArtefactualOverlap(const set<pair<string, string>>& overlapping_gene_set_pairs,
                                                                 const set<pair<string, string>>& overlapping_protein_set_pairs,
                                                                 const set<pair<string, string>>& overlapping_proteoform_set_pairs) {
   set<pair<string, string>> result;

   auto t0 = clock();
   cerr << "Comparing gene and proteoform pairs..." << endl;
   for (const auto& gene_set_pair : overlapping_gene_set_pairs) {
      set<pair<string, string>>::iterator it = overlapping_proteoform_set_pairs.find(gene_set_pair);
      if (it == overlapping_proteoform_set_pairs.end()) {
         result.insert(gene_set_pair);
      }
   }
   auto t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";

   t0 = clock();
   cerr << "Comparing protein and proteoform pairs..." << endl;
   for (const auto& protein_set_pair : overlapping_protein_set_pairs) {
      set<pair<string, string>>::iterator it = overlapping_proteoform_set_pairs.find(protein_set_pair);
      if (it == overlapping_proteoform_set_pairs.end()) {
         result.insert(protein_set_pair);
      }
   }
   t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";

   return result;
}

void reportPathwayPairsWithArtefactualOverlap(const set<pair<string, string>>& examples,
                                              const vector<string>& index_to_genes,
                                              const vector<string>& index_to_proteins,
                                              const vector<string>& index_to_proteoforms,
                                              const map<string, bitset<NUM_GENES>>& pathway_to_genes,
                                              const map<string, bitset<NUM_PROTEINS>>& pathway_to_proteins,
                                              const map<string, bitset<NUM_PROTEOFORMS>>& pathway_to_proteoforms,
                                              const set<pair<string, string>>& overlapping_gene_set_pairs,
                                              const set<pair<string, string>>& overlapping_protein_set_pairs,
                                              const set<pair<string, string>>& overlapping_proteoform_set_pairs,
                                              const string& report_file_path) {
   cerr << "Reporting pathway pairs with artefactual overlap...\n";
   ofstream report(report_file_path);

   report << "PATHWAY_1\tPATHWAY_2\t";
   report << "PATHWAY_1_GENE_SIZE\tPATHWAY_1_PROTEIN_SIZE\tPATHWAY_1_PROTEOFORM_SIZE\t";
   report << "PATHWAY_2_GENE_SIZE\tPATHWAY_2_PROTEIN_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\t";
   report << "GENE_OVERLAP\tPROTEIN_OVERLAP\tPROTEOFORM_OVERLAP\t";
   report << "OVERLAP_GENES\tOVERLAP_PROTEINS\tOVERLAP_PROTEOFORMS\n";
   for (const auto& example : examples) {
      report << example.first << " " << example.second << "\t";

      report << pathway_to_genes.at(example.first).count() << "\t" << pathway_to_proteins.at(example.first).count() << "\t" << pathway_to_proteoforms.at(example.first).count()
             << "\t";
      report << pathway_to_genes.at(example.second).count() << "\t" << pathway_to_proteins.at(example.second).count() << "\t" << pathway_to_proteoforms.at(example.second).count()
             << "\t";

      bitset<NUM_GENES> overlap_genes = pathway_to_genes.at(example.first) & pathway_to_genes.at(example.second);
      bitset<NUM_PROTEINS> overlap_proteins = pathway_to_proteins.at(example.first) & pathway_to_proteins.at(example.second);
      bitset<NUM_PROTEOFORMS> overlap_proteoforms = pathway_to_proteoforms.at(example.first) & pathway_to_proteoforms.at(example.second);

      report << overlap_genes.count() << "\t" << overlap_proteins.count() << "\t" << overlap_proteoforms.count() << "\t";
      printGeneMembers(report, overlap_genes, index_to_genes);
      printProteinMembers(report, overlap_proteins, index_to_proteins);
      printProteoformMembers(report, overlap_proteoforms, index_to_proteoforms);
      report << "\n";
   }
}

// Find pairs of modules/pathways that overlap on gene or protein network, but not in proteoform network
void doArtefactualOverlapAnalysis(const string& path_file_gene_search,
                                  const string& path_file_protein_search,
                                  const string& path_file_proteoform_search,
                                  const string& report_file_path) {
   cout << "Searching for artefactual overlap examples..." << endl;
   //    // TODO: Report pairs of disease modules
   //    // findDiseaseModulePairsWithArtefactualOverlap(path_file_gene_art_pairs, path_file_protein_art_pairs);

   vector<string> index_to_genes = loadEntities(path_file_gene_search);
   vector<string> index_to_proteins = loadEntities(path_file_protein_search);
   vector<string> index_to_proteoforms = loadEntities(path_file_proteoform_search);

   map<string, int> genes_to_index = fillMap(index_to_genes);
   map<string, int> proteins_to_index = fillMap(index_to_proteins);
   map<string, int> proteoforms_to_index = fillMap(index_to_proteoforms);

   map<string, bitset<NUM_GENES>> pathway_to_genes = loadPathwaysGeneMembers(path_file_gene_search, genes_to_index);
   map<string, bitset<NUM_PROTEINS>> pathway_to_proteins = loadPathwaysProteinMembers(path_file_protein_search, proteins_to_index);
   map<string, bitset<NUM_PROTEOFORMS>> pathway_to_proteoforms = loadPathwaysProteoformMembers(path_file_proteoform_search, proteoforms_to_index);

   // Select pathways that fullfill the requirements

   cout << "Calculating gene sets overlap..." << endl;
   set<pair<string, string>> overlapping_gene_set_pairs = findOverlappingGeneSets(pathway_to_genes, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_PATHWAY_SIZE, MAX_PATHWAY_SIZE);
   cout << "Calculating protein sets overlap..." << endl;
   set<pair<string, string>> overlapping_protein_set_pairs = findOverlappingProteinSets(pathway_to_proteins, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_PATHWAY_SIZE, MAX_PATHWAY_SIZE);
   cout << "Calculating proteoform sets..." << endl;
   set<pair<string, string>> overlapping_proteoform_set_pairs = findOverlappingProteoformSets(pathway_to_proteoforms);

   cout << "Finding examples of artifactual overlap...\n";
   set<pair<string, string>> examples = findPathwayPairsWithArtefactualOverlap(overlapping_gene_set_pairs,
                                                                               overlapping_protein_set_pairs,
                                                                               overlapping_proteoform_set_pairs);

   //TODO: Report the different states of the proteins that are shared in the protein level, but not in the proteoform level.
   reportPathwayPairsWithArtefactualOverlap(examples,
                                            index_to_genes,
                                            index_to_proteins,
                                            index_to_proteoforms,
                                            pathway_to_genes,
                                            pathway_to_proteins,
                                            pathway_to_proteoforms,
                                            overlapping_gene_set_pairs,
                                            overlapping_protein_set_pairs,
                                            overlapping_proteoform_set_pairs,
                                            report_file_path);
}