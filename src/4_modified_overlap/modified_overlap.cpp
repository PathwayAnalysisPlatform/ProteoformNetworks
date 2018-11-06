#include "modified_overlap.hpp"

using namespace std;

bool isModified(const string& proteoform) {
   smatch modification;
   return regex_search(proteoform, modification, RGX_MODIFICATION);
}

bitset<NUM_PROTEOFORMS> getSetOfModifiedProteoforms(const vector<string>& proteoforms) {
   bitset<NUM_PROTEOFORMS> modified_proteoforms;

   for (int I = 0; I < proteoforms.size(); I++) {
      if (isModified(proteoforms[I])) {
         modified_proteoforms.set(I);
      }
   }

   return modified_proteoforms;
}

void reportPathwayPairsWithKeyPTMOverlap(const string& path_file_proteoform_search, const string& report_file_path) {
   vector<string> index_to_proteoforms = loadEntities(path_file_proteoform_search);
   map<string, int> proteoforms_to_index = fillMap(index_to_proteoforms);
   map<string, bitset<NUM_PROTEOFORMS>> pathways_to_proteoforms = loadPathwaysProteoformMembers(path_file_proteoform_search, proteoforms_to_index);
   bitset<NUM_PROTEOFORMS> modified_proteoforms = getSetOfModifiedProteoforms(index_to_proteoforms);
   vector<typename map<string, bitset<NUM_PROTEOFORMS>>::const_iterator> modifiedPathways;

   // Compare all the pairs of selected pathways and select pairs that overlap only in a percentage of modified proteins
   set<pair<string, string>> examples = findOverlappingProteoformSets(pathways_to_proteoforms,
                                                                      MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE,
                                                                      MIN_PATHWAY_SIZE, MAX_PATHWAY_SIZE,
                                                                      modified_proteoforms,
                                                                      MIN_MODIFIED_ALL_MEMBERS_RATIO,
                                                                      MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

   cout << "Reporting pathway pairs with only modified overlap...\n";
   ofstream report(report_file_path);

   report << "PATHWAY_1\tPATHWAY_2\t";
   report << "PATHWAY_1_PROTEOFORM_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\tOVERLAP_SIZE\t";
   report << "OVERLAP_PROTEOFORMS\n";
   for (const auto& example : examples) {
      bitset<NUM_PROTEOFORMS> overlap_proteoforms = pathways_to_proteoforms.at(example.first) & pathways_to_proteoforms.at(example.second);

      report << example.first << " " << example.second << "\t";
      report << pathways_to_proteoforms.at(example.first).count() << "\t" << pathways_to_proteoforms.at(example.first).count() << "\t";
      report << overlap_proteoforms.count() << "\t";
      printProteoformMembers(report, overlap_proteoforms, index_to_proteoforms);

      report << "\n";
   }

   // TODO: Review the types of modifications that appear in the overlaps
}

void doModifiedOverlapAnalysis(const string& path_file_proteoform_search,
                               const string& path_file_report,
                               const string& path_file_proteins,
                               const string& path_file_proteoforms,
                               const string& path_file_modifications) {
   cout << "Searching for modified overlap examples...\n";

   // Part 1: Find examples of pathways that overlap only in modified proteins
   reportPathwayPairsWithKeyPTMOverlap(path_file_proteoform_search, path_file_report);

   // Part 2: Find examples of disease modules that overlap only in modified proteins
   //    reportDiseaseModulePairsWithKeyPTMOverlap();
}