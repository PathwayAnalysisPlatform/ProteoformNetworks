#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std;

const size_t NUM_PROTEOFORMS = 13911;

const int MIN_OVERLAP_SIZE = 1;
const int MAX_OVERLAP_SIZE = 10;

const int MIN_PATHWAY_SIZE = 1;
const int MAX_PATHWAY_SIZE = 20;

const bool SHOW_PATHWAY_NAMES = true;

const std::regex RGX_MODIFICATION{"[;,]\\d{5}"};

const float MIN_MODIFIED_ALL_MEMBERS_RATIO = 0.1;
const float MIN_MODIFIED_OVERLAP_MEMBERS_RATIO = 0.9;

vector<string> loadEntities(const string& entities_file_path) {
   ifstream entities_file(entities_file_path);
   string entity, line_leftover;
   unordered_set<string> temp_set;
   vector<string> entities;

   if (!entities_file.is_open()) {
      throw std::runtime_error("Cannot open " + entities_file_path);
   }

   while (getline(entities_file, entity, '\t')) {
      temp_set.insert(entity);
      getline(entities_file, line_leftover);
   }
   entities.assign(temp_set.begin(), temp_set.end());
   sort(entities.begin(), entities.end());
   return entities;
}

map<string, int> fillMap(const vector<string>& index_to_entities) {
   map<string, int> entities_to_index;
   for (int I = 0; I < index_to_entities.size(); I++) {
      entities_to_index.emplace(index_to_entities[I], I);
   }
   return entities_to_index;
}

map<string, bitset<NUM_PROTEOFORMS>> loadPathwaysProteoformMembers(const string& file_path, const map<string, int>& entities_to_index) {
   map<string, bitset<NUM_PROTEOFORMS>> result;
   ifstream file_search(file_path);
   string field, entity, pathway;
   bitset<NUM_PROTEOFORMS> empty_set;

   getline(file_search, field);                  // Skip csv header line
   while (getline(file_search, entity, '\t')) {  // Read the members of each pathway // Read PROTEOFORM
      getline(file_search, field, '\t');         // Read UNIPROT
      getline(file_search, field, '\t');         // Read REACTION_STID
      getline(file_search, field, '\t');         // Read REACTION_DISPLAY_NAME
      getline(file_search, pathway, '\t');       // Read PATHWAY_STID
      getline(file_search, field);               // Read PATHWAY_DISPLAY_NAME

      if (result.find(pathway) == result.end()) {
         result.emplace(pathway, empty_set);
      }
      if (entities_to_index.find(entity) != entities_to_index.end()) {
         result[pathway].set(entities_to_index.find(entity)->second);
      }
   }
   return result;
}

template <size_t set_size>
void getMembers(const bitset<set_size>& pathway_members, const vector<string>& index_to_entities) {
   for (int I = 0; I < set_size; I++) {
      if (overlap.test(I)) {
         output << index_to_entities[I];
         printed++;
      }
   }
   output << "\t";
}

bool isModified(const string& proteoform) {
   std::smatch modification;
   return std::regex_search(proteoform, modification, RGX_MODIFICATION);
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

// vector<string> getSetsWithModifiedProteoforms(const bitset<NUM_PROTEOFORMS> modified_proteoforms,
//                                               const map<string, bitset<NUM_PROTEOFORMS>> sets_to_proteoforms,
//                                               const float& ratio) {
//    vector<string> result;
//    for (const auto& set_entry : sets_to_proteoforms) {
//       float percentage = static_cast<float>((modified_proteoforms & set_entry.second).count()) / static_cast<float>(set_entry.second.count());
//       if (percentage >= ratio) {
//          result.push_back(set_entry.first);
//       }
//    }
//    return result;
// }

float getModifiedRatio(const bitset<NUM_PROTEOFORMS>& modified_proteoforms, const bitset<NUM_PROTEOFORMS> set_members) {
   return static_cast<float>((modified_proteoforms & set_members).count()) / static_cast<float>(set_members.count());
}

template <size_t set_size>
vector<typename map<string, bitset<set_size>>::const_iterator> getSetsWithModifiedProteoforms(const bitset<NUM_PROTEOFORMS>& modified_proteoforms,
                                                                                              const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_proteoforms,
                                                                                              const int& min_set_size, const int& max_set_size,
                                                                                              const float& ratio) {
   std::vector<typename map<string, bitset<set_size>>::const_iterator> result;
   for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
      int set_size = it->second.count();
      if (min_set_size <= set_size && set_size <= max_set_size) {
         float percentage = static_cast<float>((modified_proteoforms & it->second).count()) / static_cast<float>(set_size);
         if (percentage >= ratio) {
            result.push_back(it);
         }
      }
   }
}

template <size_t set_size>
set<pair<string, string>> findOverlappingPairs(const std::vector<typename map<string, bitset<set_size>>::const_iterator>& nav,
                                               const map<string, bitset<set_size>>& sets_to_members,
                                               int min_overlap, int max_overlap,
                                               int min_set_size, int max_set_size) {
   std::cerr << "elementos: " << sets_to_members.size() << ", parejas: " << (sets_to_members.size() * sets_to_members.size() - sets_to_members.size()) / 2 << "\n";
   std::cerr << "Searching for sets in the size range [" << min_set_size << ", " << max_set_size << "]\n";
   auto t0 = std::clock();
   set<pair<string, string>> result;
   for (auto vit1 = nav.begin(); vit1 != nav.end(); vit1++) {
      for (auto vit2 = std::next(vit1); vit2 != nav.end(); vit2++) {
         if (min_set_size <= (*vit1)->second.count() && (*vit1)->second.count() <= max_set_size && min_set_size <= (*vit2)->second.count() && (*vit2)->second.count() <= max_set_size) {
            bitset<set_size> overlap = (*vit1)->second & (*vit2)->second;
            if (min_overlap <= overlap.count() && overlap.count() <= max_overlap) {
               result.emplace((*vit1)->first, (*vit2)->first);
            }
         }
      }
   }
   auto t1 = std::clock();
   std::cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
   return result;
}

set<pair<string, string>> findKeyPTMOverlappingPairs(const set<pair<string, string>>& overlaping_proteoform_set_pairs,
                                                     const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_proteoforms,
                                                     const bitset<NUM_PROTEOFORMS>& modified_proteoforms,
                                                     const float& min_modified_overlap_members_ratio) {
   set<pair<string, string>> result;

   auto t0 = std::clock();
   std::cerr << "Comparing protein and proteoform pairs..." << endl;
   for (const auto& set_pair : overlaping_proteoform_set_pairs) {
      bitset<NUM_PROTEOFORMS> overlap = sets_to_proteoforms.at(set_pair.first) & sets_to_proteoforms.at(set_pair.second);
      float percentage = static_cast<float>((modified_proteoforms & overlap).count()) / static_cast<float>(overlap.count());
      if (percentage >= min_modified_overlap_members_ratio) {
         result.insert(set_pair);
      }
   }
   auto t1 = std::clock();
   std::cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";

   return result;
}

template <size_t set_size>
void printMembers(std::ofstream& output, const bitset<set_size>& overlap, const vector<string>& index_to_entities) {
   int printed = 0;
   int total = overlap.count();
   for (int I = 0; I < set_size; I++) {
      if (overlap.test(I)) {
         output << index_to_entities[I];
         printed++;
         if (printed != total) {
            output << ";";
         }
      }
   }
   output << "\t";
}

void reportPathwayPairsWithKeyPTMOverlap(const set<pair<string, string>>& examples,
                                         const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_proteoforms,
                                         const vector<string>& index_to_proteoforms,
                                         const string& report_file_path) {
   std::cerr << "Reporting pathway pairs with only modified overlap...\n";
   ofstream report(report_file_path);

   report << "PATHWAY_1\tPATHWAY_2\t";
   report << "PATHWAY_1_PROTEOFORM_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\tOVERLAP_SIZE\t";
   report << "OVERLAP_PROTEOFORMS\n";
   for (const auto& example : examples) {
      bitset<NUM_PROTEOFORMS> overlap_proteoforms = sets_to_proteoforms.at(example.first) & sets_to_proteoforms.at(example.second);

      report << example.first << " " << example.second << "\t";
      report << sets_to_proteoforms.at(example.first).count() << "\t" << sets_to_proteoforms.at(example.first).count() << "\t";
      report << overlap_proteoforms.count() << "\t";
      printMembers(report, overlap_proteoforms, index_to_proteoforms);

      report << "\n";
   }
}

void doKeyPTMOverlapAnalysis(const string& path_file_gene_search, const string& path_file_protein_search, const string& path_file_proteoform_search,
                             const string& report_file_path) {
   vector<string> index_to_proteoforms = loadEntities(path_file_proteoform_search);
   map<string, int> proteoforms_to_index = fillMap(index_to_proteoforms);
   map<string, bitset<NUM_PROTEOFORMS>> pathways_to_proteoforms = loadPathwaysProteoformMembers(path_file_proteoform_search, proteoforms_to_index);
   bitset<NUM_PROTEOFORMS> modified_proteoforms = getSetOfModifiedProteoforms(index_to_proteoforms);
   vector<typename map<string, bitset<NUM_PROTEOFORMS>>::const_iterator> modifiedPathways;

   modifiedPathways = getSetsWithModifiedProteoforms(modified_proteoforms,
                                                     pathways_to_proteoforms,
                                                     MIN_PATHWAY_SIZE, MAX_PATHWAY_SIZE,
                                                     MIN_MODIFIED_ALL_MEMBERS_RATIO);

   // Compare all the pairs of selected pathways and select pairs that overlap only in a percentage of modified proteins
   set<pair<string, string>> overlaping_proteoform_set_pairs = findOverlappingPairs(modifiedPathways,
                                                                                    pathways_to_proteoforms,
                                                                                    MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE,
                                                                                    MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

   set<pair<string, string>> examples = findKeyPTMOverlappingPairs(overlaping_proteoform_set_pairs,
                                                                   pathways_to_proteoforms,
                                                                   modified_proteoforms, 
                                                                   MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

   reportPathwayPairsWithKeyPTMOverlap(examples, pathways_to_proteoforms, index_to_proteoforms, report_file_path);
}