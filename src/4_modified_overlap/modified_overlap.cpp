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
                        const unordered_map<string, string>& pathways_to_names,
                        const unordered_map<string, bitset<NUM_PROTEOFORMS>>& pathways_to_proteoforms,
                        const vector<string> index_to_proteoforms) {
   for (const auto& example : examples) {
      output << example.first.first << "\t" << example.first.second << "\t" << pathways_to_names.at(example.first.first) << "\t"
             << pathways_to_names.at(example.first.second) << "\t";
      output << pathways_to_proteoforms.at(example.first.first).count() << "\t" << pathways_to_proteoforms.at(example.first.second).count() << "\t";
      output << example.second.count() << "\t";
      printMembers(output, example.second, index_to_proteoforms);
      output << "\n";
   }
}

void writePathwayReport(ofstream& output,
                        const map<pair<string, string>, bitset<NUM_PROTEOFORMS>>& examples,
                        const unordered_map<string, string>& pathways_to_names,
                        const unordered_map<string, bitset<NUM_PROTEOFORMS>>& pathways_to_proteoforms,
                        const vector<string> index_to_proteoforms) {
   output << "PATHWAY_1\tPATHWAY_2\tPATHWAY_1_NAME\tPATHWAY_2_NAME\t";
   output << "PATHWAY_1_PROTEOFORM_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\tOVERLAP_SIZE\t";
   output << "OVERLAP_PROTEOFORMS\n";
   writeReportRecords(output, examples, pathways_to_names, pathways_to_proteoforms, index_to_proteoforms);
}

void writePhenotypeReport(ofstream& output,
                          const map<pair<string, string>, bitset<NUM_PROTEOFORMS>>& examples,
                          const unordered_map<string, string>& pathways_to_names,
                          const unordered_map<string, bitset<NUM_PROTEOFORMS>>& pathways_to_proteoforms,
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

void writeFrequencies(const string& modifications_file_path, const string& proteins_file_path, const string& proteoforms_file_path,
                      const map<string, int>& mod_to_freq, const map<string, int>& proteins_to_freq, const map<string, int>& proteoforms_to_freq) {
   ofstream modifications_file(modifications_file_path);
   ofstream proteins_file(proteins_file_path);
   ofstream proteoforms_file(proteoforms_file_path);

   modifications_file << "MODIFICATION\tFREQUENCY\n";
   for (const auto& modification : mod_to_freq) {
      modifications_file << modification.first << "\t" << modification.second << "\n";
   }
   cerr << "Finished writing modification frequencies.\n";

   proteins_file << "PROTEIN\tFREQUENCY\n";
   for (const auto& protein : proteins_to_freq) {
      proteins_file << protein.first << "\t" << protein.second << "\n";
   }
   cerr << "Finished writing protein frequencies.\n";

   proteoforms_file << "PROTEOFORM\tFREQUENCY\n";
   for (const auto& proteoform : proteoforms_to_freq) {
      proteoforms_file << proteoform.first << "\t" << proteoform.second << "\n";
   }
   cerr << "Finished writing proteoform frequencies.\n";
}

Frequencies getFrequencies(const map<pair<string, string>, bitset<NUM_PROTEOFORMS>>& proteoform_overlap_pairs,
                           const vector<string>& index_to_proteoforms) {
   Frequencies frequencies;
   for (const auto& overlap_pair : proteoform_overlap_pairs) {
      // For each overlap set
      for (int I = 0; I < overlap_pair.second.size(); I++) {
         if (overlap_pair.second.test(I)) {
            string accession = getAccession(index_to_proteoforms[I]);
            if (frequencies.proteins.find(accession) == frequencies.proteins.end()) {
               frequencies.proteins.emplace(accession, 0);
            }
            frequencies.proteins[accession]++;

            if (frequencies.proteoforms.find(index_to_proteoforms[I]) == frequencies.proteoforms.end()) {
               frequencies.proteoforms.emplace(index_to_proteoforms[I], 0);
            }
            frequencies.proteoforms[index_to_proteoforms[I]]++;

            for (const auto& modification : getModifications(index_to_proteoforms[I])) {
               if (frequencies.modifications.find(modification) == frequencies.modifications.end()) {
                  frequencies.modifications[modification] = 0;
               }
               frequencies.modifications[modification]++;
            }
         }
      }
   }
   return frequencies;
}

void plotFrequencies(string report_file_path, string modifications_file_path, string proteins_file_path, string proteoforms_file_path) {
   string command = "Rscript src/4_modified_overlap/modified_overlap.R " + report_file_path + " " +
                    modifications_file_path + " " + proteins_file_path + " " + proteoforms_file_path;
   cerr << "Plotting modified overlap frequencies: " << command << endl;
   system(command.c_str());
}

void reportPathwayPairs(const string& path_file_proteoform_search,
                        const string& report_file_path,
                        const string& modifications_file_path,
                        const string& proteins_file_path,
                        const string& proteoforms_file_path) {
   const auto [index_to_proteoforms, proteoforms_to_index] = loadEntities(path_file_proteoform_search);
   const auto pathways_to_names = loadPathwayNames(path_file_proteoform_search);
   const auto pathways_to_proteoforms = loadProteoformSets(path_file_proteoform_search, proteoforms_to_index, true);
   const bitset<NUM_PROTEOFORMS> modified_proteoforms = getSetOfModifiedProteoforms(index_to_proteoforms);

   cout << "Reporting pathway pairs with only modified overlap...\n";

   // Compare all the pairs of selected pathways and select pairs that overlap only in a percentage of modified proteins
   const auto& examples =
       findOverlappingProteoformSets(pathways_to_proteoforms, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE,
                                     modified_proteoforms, MIN_MODIFIED_ALL_MEMBERS_RATIO, MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

   ofstream report(report_file_path);
   writePathwayReport(report, examples, pathways_to_names, pathways_to_proteoforms, index_to_proteoforms);

   const auto [mod_to_freq, proteins_to_freq, proteoforms_to_freq] = getFrequencies(examples, index_to_proteoforms);
   writeFrequencies(modifications_file_path, proteins_file_path, proteoforms_file_path, mod_to_freq, proteins_to_freq, proteoforms_to_freq);
   plotFrequencies(report_file_path, modifications_file_path, proteins_file_path, proteoforms_file_path);
}

void reportPhenotypePairs(const string& path_file_proteoform_search,
                          const string& path_file_PheGenI_full,
                          const std::string& path_file_mapping_proteins_to_genes,
                          const string& path_file_report_trait,
                          const string& path_file_modified_overlap_trait_proteins,
                          const string& path_file_modified_overlap_trait_proteoforms,
                          const string& path_file_modified_overlap_trait_modifications) {
   const auto [index_to_proteoforms, proteoforms_to_index] = loadEntities(path_file_proteoform_search);

   // TODO: Load trait gene sets

   // TODO: Convert gene sets to proteoform sets
}

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
                const std::string& path_file_modified_overlap_trait_modifications) {
   cout << "Searching for modified overlap examples...\n";

   // Part 1: Find examples of pathways that overlap only in modified proteins
   reportPathwayPairs(path_file_proteoform_search,
                      path_file_report_pathway,
                      path_file_modified_overlap_pathway_modifications,
                      path_file_modified_overlap_pathway_proteins,
                      path_file_modified_overlap_pathway_proteoforms);

   // Part 2: Find examples of disease modules that overlap only in modified proteins
   reportPhenotypePairs(path_file_proteoform_search,
                        path_file_PheGenI_full,
                        path_file_mapping_proteins_to_genes,
                        path_file_report_trait,
                        path_file_modified_overlap_trait_proteins,
                        path_file_modified_overlap_trait_proteoforms,
                        path_file_modified_overlap_trait_modifications);
}

}  // namespace modified_overlap
