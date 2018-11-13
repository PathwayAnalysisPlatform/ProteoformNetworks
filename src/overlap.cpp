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

map<string, string> loadPathwayNames(const string& path_search_file) {
   map<string, string> result;
   ifstream file_search(path_search_file);
   string field, pathway, pathway_name;
   int field_cont = 0;
   int pathway_stid_column = 0;

   while (getline(file_search, field, '\t')) {
      if (field == "PATHWAY_STID") {
         pathway_stid_column = field_cont;
         break;
      }
      field_cont++;
   }
   getline(file_search, field);  // Skip header line leftoever

   while (getline(file_search, field, '\t')) {                        // While there are records in the file
      for (int column = 1; column < pathway_stid_column; column++) {  // Skip fields before the pathway_stid
         getline(file_search, field, '\t');
      }
      getline(file_search, pathway, '\t');  // Read PATHWAY_STID
      getline(file_search, pathway_name);   // Read PATHWAY_DISPLAY_NAME
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
map<pair<string, string>, bitset<total_num_entities>> findOverlappingPairs(const map<string, bitset<total_num_entities>>& sets_to_members,
                                                                           const int& min_overlap, const int& max_overlap,
                                                                           const int& min_set_size, const int& max_set_size) {
   map<pair<string, string>, bitset<total_num_entities>> result;
   vector<typename map<string, bitset<total_num_entities>>::const_iterator> nav;

   for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
      int set_size = it->second.count();
      if (min_set_size <= set_size && set_size <= max_set_size) {
         nav.push_back(it);
      }
   }

   cerr << "elementos: " << nav.size() << ", parejas: " << (nav.size() * nav.size() - nav.size()) / 2 << "\n";
   auto t0 = clock();

   for (auto vit1 = nav.begin(); vit1 != nav.end(); vit1++) {
      for (auto vit2 = next(vit1); vit2 != nav.end(); vit2++) {
         bitset<total_num_entities> overlap = (*vit1)->second & (*vit2)->second;
         if (min_overlap <= overlap.count() && overlap.count() <= max_overlap) {
            result.emplace(make_pair((*vit1)->first, (*vit2)->first), overlap);
         }
      }
   }
   auto t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
   return result;
}

map<pair<string, string>, bitset<NUM_GENES>> findOverlappingGeneSets(const map<string, bitset<NUM_GENES>>& sets_to_members,
                                                                     const int& min_overlap, const int& max_overlap,
                                                                     const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members, min_overlap, max_overlap, min_set_size, max_set_size);
}

map<pair<string, string>, bitset<NUM_PROTEINS>> findOverlappingProteinSets(const map<string, bitset<NUM_PROTEINS>>& sets_to_members,
                                                                           const int& min_overlap, const int& max_overlap,
                                                                           const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members,
                               min_overlap, max_overlap,
                               min_set_size, max_set_size);
}

map<pair<string, string>, bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                                 const int& min_overlap, const int& max_overlap,
                                                                                 const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members,
                               min_overlap, max_overlap,
                               min_set_size, max_set_size);
}

//Version without set size limits or overlap size limits
map<pair<string, string>, bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members) {
   return findOverlappingPairs(sets_to_members,
                               1, NUM_PROTEOFORMS,
                               1, NUM_PROTEOFORMS);
}

// This version includes modified ratio of the proteoform members of the set, and another ratio for the overlapping proteoforms.
map<pair<string, string>, bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                                 const int& min_overlap, const int& max_overlap,
                                                                                 const int& min_set_size, const int& max_set_size,
                                                                                 const bitset<NUM_PROTEOFORMS>& modified_proteoforms,
                                                                                 const float& min_all_modified_ratio,
                                                                                 const float& min_overlap_modified_ratio) {
   vector<typename map<string, bitset<NUM_PROTEOFORMS>>::const_iterator> nav;
   for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
      int set_size = it->second.count();
      if (min_set_size <= set_size && set_size <= max_set_size) {
        //  cerr << "MIN_SET_SIZE: " << min_set_size << " SET_SIZE: " << set_size << " MAX_SET_SIZE: " << max_set_size << "\n";
         float percentage = static_cast<float>((modified_proteoforms & it->second).count()) / static_cast<float>(set_size);
         if (percentage >= min_all_modified_ratio) {
            nav.push_back(it);
         }
      }
   }

   map<pair<string, string>, bitset<NUM_PROTEOFORMS>> result;
   cerr << "elementos: " << nav.size() << ", parejas: " << (nav.size() * nav.size() - nav.size()) / 2 << "\n";
   auto t0 = clock();
   for (auto vit1 = nav.begin(); vit1 != nav.end(); vit1++) {
      for (auto vit2 = next(vit1); vit2 != nav.end(); vit2++) {
         bitset<NUM_PROTEOFORMS> overlap = (*vit1)->second & (*vit2)->second;
         if (min_overlap <= overlap.count() && overlap.count() <= max_overlap) {
            float percentage = static_cast<float>((modified_proteoforms & overlap).count()) / static_cast<float>(overlap.count());
            if (percentage >= min_overlap_modified_ratio) {
               result.emplace(make_pair((*vit1)->first, (*vit2)->first), overlap);
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
   output << "[";
   for (int I = 0; I < total_num_entities; I++) {
      if (entity_set.test(I)) {
         output << "\"" << index_to_entities[I] << "\"";
         printed++;
         if (printed != total) {
            output << ",";
         }
      }
   }
   output << "]";
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

template <size_t total_num_entities>
set<string> getEntityStrings(const bitset<total_num_entities>& entity_set, const vector<string>& index_to_entities) {
   set<string> result;
   for (int I = 0; I < total_num_entities; I++) {
      if (entity_set.test(I)) {
         result.insert(index_to_entities[I]);
      }
   }
   return result;
}

set<string> getGeneStrings(const bitset<NUM_GENES>& gene_set, const vector<string>& index_to_genes) {
   return getEntityStrings(gene_set, index_to_genes);
}

set<string> getProteinStrings(const bitset<NUM_PROTEINS>& protein_set, const vector<string>& index_to_proteins) {
   return getEntityStrings(protein_set, index_to_proteins);
}

set<string> getProteoformStrings(const bitset<NUM_PROTEOFORMS>& proteoform_set, const vector<string>& index_to_proteoforms) {
   return getEntityStrings(proteoform_set, index_to_proteoforms);
}


string getAccession(string proteoform) {
   smatch match_end_of_accession;
   if (!regex_search(proteoform, match_end_of_accession, RGX_ACCESSION_DELIMITER)) {
      return proteoform;
   }
   return proteoform.substr(0, match_end_of_accession.position(0));
}