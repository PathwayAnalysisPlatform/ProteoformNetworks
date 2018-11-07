#include "overlap.hpp"

using namespace std;

vector<string> getIndexToEntities(const string& entities_file_path) {
   ifstream entities_file(entities_file_path);
   string entity, line_leftover;
   unordered_set<string> temp_set;
   vector<string> index_to_entities;

   if (!entities_file.is_open()) {
      throw runtime_error("Cannot open " + entities_file_path);
   }
   
   while (getline(entities_file, entity, '\t')) {
      temp_set.insert(entity);
      getline(entities_file, line_leftover);
   }
   index_to_entities.assign(temp_set.begin(), temp_set.end());
   sort(index_to_entities.begin(), index_to_entities.end());

   return index_to_entities;
}

map<string, int> getEntitiesToIndex(const vector<string>& index_to_entities) {
   map<string, int> entities_to_index;
   for (int I = 0; I < index_to_entities.size(); I++) {
      entities_to_index.emplace(index_to_entities[I], I);
   }
   return entities_to_index;
}

Entities_bimap loadEntities(const string& entities_file_path) {
    auto index_to_entities = getIndexToEntities(entities_file_path);
    auto entities_to_index = getEntitiesToIndex(index_to_entities);
   return {index_to_entities, entities_to_index};
}

map<string, string> loadPathwayNames(const string& path_protein_search_file) {
   map<string, string> result;
   ifstream file_search(path_protein_search_file);
   string field, pathway, pathway_name;

   getline(file_search, field);                 // Skip csv header line
   while (getline(file_search, field, '\t')) {  // Read UNIPROT
      getline(file_search, field, '\t');        // Read REACTION_STID
      getline(file_search, field, '\t');        // Read REACTION_DISPLAY_NAME
      getline(file_search, pathway, '\t');      // Read PATHWAY_STID
      getline(file_search, pathway_name);       // Read PATHWAY_DISPLAY_NAME

      result[pathway] = pathway_name;
   }
   return result;
}

map<string, bitset<NUM_GENES>> loadPathwaysGeneMembers(const string& file_path, const map<string, int>& entities_to_index) {
   map<string, bitset<NUM_GENES>> result;
   ifstream file_search(file_path);
   string field, gene, pathway;
   bitset<NUM_GENES> empty_set;

   getline(file_search, field);                // Skip csv header line
   while (getline(file_search, gene, '\t')) {  // Read the members of each pathway // Read GENE
      getline(file_search, field, '\t');       // Read UNIPROT
      getline(file_search, field, '\t');       // Read REACTION_STID
      getline(file_search, field, '\t');       // Read REACTION_DISPLAY_NAME
      getline(file_search, pathway, '\t');     // Read PATHWAY_STID
      getline(file_search, field);             // Read PATHWAY_DISPLAY_NAME

      if (result.find(pathway) == result.end()) {
         result.emplace(pathway, empty_set);
      }
      if (entities_to_index.find(gene) != entities_to_index.end()) {
         result[pathway].set(entities_to_index.find(gene)->second);
      }
   }
   return result;
}

map<string, bitset<NUM_PROTEINS>> loadPathwaysProteinMembers(const string& file_path, const map<string, int>& entities_to_index) {
   map<string, bitset<NUM_PROTEINS>> result;
   ifstream file_search(file_path);
   string field, entity, pathway;
   bitset<NUM_PROTEINS> empty_set;

   getline(file_search, field);                  // Skip csv header line
   while (getline(file_search, entity, '\t')) {  // Read the members of each pathway // Read UNIPROT
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

// Version with set size and overlap size limits
template <size_t total_num_entities>
set<pair<string, string>> findOverlappingPairs(const map<string, bitset<total_num_entities>>& sets_to_members,
                                               const int& min_overlap, const int& max_overlap,
                                               const int& min_set_size, const int& max_set_size) {
   vector<typename map<string, bitset<total_num_entities>>::const_iterator> nav;
   for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
      int set_size = it->second.count();
      if (min_set_size <= set_size && set_size <= max_set_size) {
         nav.push_back(it);
      }
   }

   cerr << "elementos: " << nav.size() << ", parejas: " << (nav.size() * nav.size() - nav.size()) / 2 << "\n";
   auto t0 = clock();
   set<pair<string, string>> result;
   for (auto vit1 = nav.begin(); vit1 != nav.end(); vit1++) {
      for (auto vit2 = next(vit1); vit2 != nav.end(); vit2++) {
         bitset<total_num_entities> overlap = (*vit1)->second & (*vit2)->second;
         if (min_overlap <= overlap.count() && overlap.count() <= max_overlap) {
            result.emplace((*vit1)->first, (*vit2)->first);
         }
      }
   }
   auto t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
   return result;
}

set<pair<string, string>> findOverlappingGeneSets(const map<string, bitset<NUM_GENES>>& sets_to_members,
                                                  const int& min_overlap, const int& max_overlap,
                                                  const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members,
                               min_overlap, max_overlap,
                               min_set_size, max_set_size);
}

set<pair<string, string>> findOverlappingProteinSets(const map<string, bitset<NUM_PROTEINS>>& sets_to_members,
                                                     const int& min_overlap, const int& max_overlap,
                                                     const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members,
                               min_overlap, max_overlap,
                               min_set_size, max_set_size);
}

set<pair<string, string>> findOverlappingProteoformSets(const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                        const int& min_overlap, const int& max_overlap,
                                                        const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members,
                               min_overlap, max_overlap,
                               min_set_size, max_set_size);
}

//Version without set size limits or overlap size limits
set<pair<string, string>> findOverlappingProteoformSets(const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members) {
   return findOverlappingPairs(sets_to_members,
                               1, NUM_PROTEOFORMS,
                               1, NUM_PROTEOFORMS);
}

// This version includes modified ratio of the proteoform members of the set, and another ratio for the overlapping proteoforms.
set<pair<string, string>> findOverlappingProteoformSets(const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                        const int& min_overlap, const int& max_overlap,
                                                        const int& min_set_size, const int& max_set_size,
                                                        const bitset<NUM_PROTEOFORMS>& modified_proteoforms,
                                                        const float& min_all_modified_ratio,
                                                        const float& min_overlap_modified_ratio) {
   vector<typename map<string, bitset<NUM_PROTEOFORMS>>::const_iterator> nav;
   for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
      int set_size = it->second.count();
      if (min_set_size <= set_size && set_size <= max_set_size) {
         float percentage = static_cast<float>((modified_proteoforms & it->second).count()) / static_cast<float>(set_size);
         if (percentage >= min_all_modified_ratio) {
            nav.push_back(it);
         }
      }
   }

   cerr << "elementos: " << nav.size() << ", parejas: " << (nav.size() * nav.size() - nav.size()) / 2 << "\n";
   auto t0 = clock();
   set<pair<string, string>> result;
   for (auto vit1 = nav.begin(); vit1 != nav.end(); vit1++) {
      for (auto vit2 = next(vit1); vit2 != nav.end(); vit2++) {
         bitset<NUM_PROTEOFORMS> overlap = (*vit1)->second & (*vit2)->second;
         if (min_overlap <= overlap.count() && overlap.count() <= max_overlap) {
            float percentage = static_cast<float>((modified_proteoforms & overlap).count()) / static_cast<float>(overlap.count());
            if (percentage >= min_overlap_modified_ratio) {
               result.emplace((*vit1)->first, (*vit2)->first);
            }
         }
      }
   }
   auto t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
   return result;
}

template <size_t total_num_entities>
void printMembers(ofstream& output, const bitset<total_num_entities>& entity_set, const vector<string>& index_to_entities) {
   int printed = 0;
   int total = entity_set.count();
   for (int I = 0; I < total_num_entities; I++) {
      if (entity_set.test(I)) {
         output << index_to_entities[I];
         printed++;
         if (printed != total) {
            output << ";";
         }
      }
   }
   output << "\t";
}

void printGeneMembers(ofstream& output, const bitset<NUM_GENES>& gene_set, const vector<string>& index_to_entities) {
   printMembers(output, gene_set, index_to_entities);
}
void printProteinMembers(ofstream& output, const bitset<NUM_PROTEINS>& protein_set, const vector<string>& index_to_entities) {
   printMembers(output, protein_set, index_to_entities);
}
void printProteoformMembers(ofstream& output, const bitset<NUM_PROTEOFORMS>& proteoform_set, const vector<string>& index_to_entities) {
   printMembers(output, proteoform_set, index_to_entities);
}