#include "modified_overlap.hpp"

using namespace std;

namespace modified_overlap {

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

void writeReportRecords(ofstream& output,
                        const map<pair<string, string>, bitset<NUM_PROTEOFORMS>>& examples,
                        const map<string, string>& pathways_to_names,
                        const map<string, bitset<NUM_PROTEOFORMS>>& pathways_to_proteoforms,
                        const vector<string> index_to_proteoforms) {
   for (const auto& example : examples) {
      output << example.first.first << "\t" << example.first.second << "\t" << pathways_to_names.at(example.first.first) << "\t" << pathways_to_names.at(example.first.second) << "\t";
      output << pathways_to_proteoforms.at(example.first.first).count() << "\t" << pathways_to_proteoforms.at(example.first.second).count() << "\t";
      output << example.second.count() << "\t";
      printProteoformMembers(output, example.second, index_to_proteoforms);
      output << "\n";
   }
}

void writePathwayReport(ofstream& output,
                        const map<pair<string, string>, bitset<NUM_PROTEOFORMS>>& examples,
                        const map<string, string>& pathways_to_names,
                        const map<string, bitset<NUM_PROTEOFORMS>>& pathways_to_proteoforms,
                        const vector<string> index_to_proteoforms) {
   output << "PATHWAY_1\tPATHWAY_2\tPATHWAY_1_NAME\tPATHWAY_2_NAME\t";
   output << "PATHWAY_1_PROTEOFORM_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\tOVERLAP_SIZE\t";
   output << "OVERLAP_PROTEOFORMS\n";
   writeReportRecords(output, examples, pathways_to_names, pathways_to_proteoforms, index_to_proteoforms);
}

void writePhenotypeReport(ofstream& output,
                          const map<pair<string, string>, bitset<NUM_PROTEOFORMS>>& examples,
                          const map<string, string>& pathways_to_names,
                          const map<string, bitset<NUM_PROTEOFORMS>>& pathways_to_proteoforms,
                          const vector<string> index_to_proteoforms) {
   output << "PHENOTYPE_1\tPHENOTYPE_1_2\tPHENOTYPE_1_NAME\tPHENOTYPE_2_NAME\t";
   output << "PHENOTYPE_1_PROTEOFORM_SIZE\tPHENOTYPE_2_PROTEOFORM_SIZE\tOVERLAP_SIZE\t";
   output << "OVERLAP_PROTEOFORMS\n";
   writeReportRecords(output, examples, pathways_to_names, pathways_to_proteoforms, index_to_proteoforms);
}

set<string> getModifications(string proteoform) {
   set<string> modifications;
   sregex_iterator it(proteoform.begin(), proteoform.end(), RGX_MODIFICATION);
   sregex_iterator end;
   while (it != end) {
      if (it->str().find(';') || it->str().find(',')) {
         modifications.insert(it->str().substr(1));
      } else {
         modifications.insert((*it)[0]);
      }
      it++;
   }
   return modifications;
}

map<string, int> getModificationsFrequencies(const map<pair<string, string>, bitset<NUM_PROTEOFORMS>>& proteoform_overlap_pairs, const vector<string>& index_to_proteoforms) {
   map<string, int> mod_to_freq;
   for (const auto& overlap_pair : proteoform_overlap_pairs) {
      // For each overlap set
      for (int I = 0; I < overlap_pair.second.size(); I++) {
         if (overlap_pair.second.test(I)) {
            for (const auto& modification : getModifications(index_to_proteoforms[I])) {
               if (mod_to_freq.find(modification) == mod_to_freq.end()) {
                  mod_to_freq[modification] = 0;
               }
               mod_to_freq[modification]++;
            }
         }
      }
   }
   return mod_to_freq;
}

void reportPathwayPairs(const string& path_file_proteoform_search, const string& report_file_path, const string& modification_file_path) {
   const auto [index_to_proteoforms, proteoforms_to_index] = loadEntities(path_file_proteoform_search);
   const map<string, string> pathways_to_names = loadPathwayNames(path_file_proteoform_search);
   const map<string, bitset<NUM_PROTEOFORMS>> pathways_to_proteoforms = loadPathwaysProteoformMembers(path_file_proteoform_search, proteoforms_to_index);
   const bitset<NUM_PROTEOFORMS> modified_proteoforms = getSetOfModifiedProteoforms(index_to_proteoforms);

   cout << "Reporting pathway pairs with only modified overlap...\n";

   // Compare all the pairs of selected pathways and select pairs that overlap only in a percentage of modified proteins
   const auto& examples = findOverlappingProteoformSets(pathways_to_proteoforms,
                                                        MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE,
                                                        MIN_PATHWAY_SIZE, MAX_PATHWAY_SIZE,
                                                        modified_proteoforms,
                                                        MIN_MODIFIED_ALL_MEMBERS_RATIO,
                                                        MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

   ofstream report(report_file_path);
   writePathwayReport(report, examples, pathways_to_names, pathways_to_proteoforms, index_to_proteoforms);

   const auto mod_to_freq = getModificationsFrequencies(examples, index_to_proteoforms);
   ofstream modification_file(modification_file_path);
   modification_file << "MODIFICATION\tFREQUENCY\n";
   for (const auto& modification : mod_to_freq) {
       modification_file << modification.first << "\t" << modification.second << "\n";
   }
}

void reportPhenotypePairs() {
}

void doAnalysis(const string& path_file_proteoform_search,
                const string& path_file_report,
                const string& path_file_proteins,
                const string& path_file_proteoforms,
                const string& path_file_modifications) {
   cout << "Searching for modified overlap examples...\n";

   // Part 1: Find examples of pathways that overlap only in modified proteins
   reportPathwayPairs(path_file_proteoform_search, path_file_report, path_file_modifications);

   // Part 2: Find examples of disease modules that overlap only in modified proteins
   reportPhenotypePairs();
}

}  // namespace modified_overlap
