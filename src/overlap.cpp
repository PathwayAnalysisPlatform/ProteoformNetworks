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
   index_to_entities = convert(temp_set);

   return index_to_entities;
}

unordered_map<string, int> getEntitiesToIndex(const vector<string>& index_to_entities) {
   unordered_map<string, int> entities_to_index;
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

unordered_map<string, string> loadPathwayNames(const string& path_search_file) {
   unordered_map<string, string> result;
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

unordered_map<string, bitset<NUM_GENES>> loadGeneSets(const string& file_path, const unordered_map<string, int>& entities_to_index, bool pathways) {
   unordered_map<string, bitset<NUM_GENES>> result;
   ifstream file_search(file_path);
   string field, gene, reaction, pathway;

   getline(file_search, field);                // Skip csv header line
   while (getline(file_search, gene, '\t')) {  // Read the members of each gene_set // Read GENE
      getline(file_search, field, '\t');       // Read UNIPROT
      getline(file_search, reaction, '\t');    // Read REACTION_STID
      getline(file_search, field, '\t');       // Read REACTION_DISPLAY_NAME
      getline(file_search, pathway, '\t');     // Read PATHWAY_STID
      getline(file_search, field);             // Read PATHWAY_DISPLAY_NAME

      result[pathways ? pathway : reaction].set(entities_to_index.find(gene)->second);
   }
   return result;
}

unordered_map<string, bitset<NUM_PROTEINS>> loadProteinSets(const string& file_path, const unordered_map<string, int>& entities_to_index, bool pathways) {
   unordered_map<string, bitset<NUM_PROTEINS>> result;
   ifstream file_search(file_path);
   string field, entity, reaction, pathway;

   getline(file_search, field);                  // Skip csv header line
   while (getline(file_search, entity, '\t')) {  // Read the members of each gene_set // Read UNIPROT
      getline(file_search, reaction, '\t');      // Read REACTION_STID
      getline(file_search, field, '\t');         // Read REACTION_DISPLAY_NAME
      getline(file_search, pathway, '\t');       // Read PATHWAY_STID
      getline(file_search, field);               // Read PATHWAY_DISPLAY_NAME

      result[pathways ? pathway : reaction].set(entities_to_index.find(entity)->second);
   }
   return result;
}

unordered_map<string, bitset<NUM_PROTEOFORMS>> loadProteoformSets(const string& file_path, const unordered_map<string, int>& entities_to_index, bool pathways) {
   unordered_map<string, bitset<NUM_PROTEOFORMS>> result;
   ifstream file_search(file_path);
   string field, entity, reaction, pathway;

   getline(file_search, field);                  // Skip csv header line
   while (getline(file_search, entity, '\t')) {  // Read the members of each gene_set // Read PROTEOFORM
      getline(file_search, field, '\t');         // Read UNIPROT
      getline(file_search, reaction, '\t');      // Read REACTION_STID
      getline(file_search, field, '\t');         // Read REACTION_DISPLAY_NAME
      getline(file_search, pathway, '\t');       // Read PATHWAY_STID
      getline(file_search, field);               // Read PATHWAY_DISPLAY_NAME

      result[pathways ? pathway : reaction].set(entities_to_index.find(entity)->second);
   }
   return result;
}

// Version with set size and overlap size limits
template <size_t total_num_entities>
map<pair<string, string>, bitset<total_num_entities>> findOverlappingPairs(const unordered_map<string, bitset<total_num_entities>>& sets_to_members,
                                                                           const int& min_overlap, const int& max_overlap,
                                                                           const int& min_set_size, const int& max_set_size) {
   map<pair<string, string>, bitset<total_num_entities>> result;
   vector<typename unordered_map<string, bitset<total_num_entities>>::const_iterator> nav;

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

map<pair<string, string>, bitset<NUM_GENES>> findOverlappingGeneSets(const unordered_map<string, bitset<NUM_GENES>>& sets_to_members,
                                                                     const int& min_overlap, const int& max_overlap,
                                                                     const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members, min_overlap, max_overlap, min_set_size, max_set_size);
}

map<pair<string, string>, bitset<NUM_PROTEINS>> findOverlappingProteinSets(const unordered_map<string, bitset<NUM_PROTEINS>>& sets_to_members,
                                                                           const int& min_overlap, const int& max_overlap,
                                                                           const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members,
                               min_overlap, max_overlap,
                               min_set_size, max_set_size);
}

map<pair<string, string>, bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const unordered_map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                                 const int& min_overlap, const int& max_overlap,
                                                                                 const int& min_set_size, const int& max_set_size) {
   return findOverlappingPairs(sets_to_members,
                               min_overlap, max_overlap,
                               min_set_size, max_set_size);
}

//Version without set size limits or overlap size limits
map<pair<string, string>, bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const unordered_map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members) {
   return findOverlappingPairs(sets_to_members,
                               1, NUM_PROTEOFORMS,
                               1, NUM_PROTEOFORMS);
}

// This version includes modified ratio of the proteoform members of the set, and another ratio for the overlapping proteoforms.
map<pair<string, string>, bitset<NUM_PROTEOFORMS>> findOverlappingProteoformSets(const unordered_map<string, bitset<NUM_PROTEOFORMS>>& sets_to_members,
                                                                                 const int& min_overlap, const int& max_overlap,
                                                                                 const int& min_set_size, const int& max_set_size,
                                                                                 const bitset<NUM_PROTEOFORMS>& modified_proteoforms,
                                                                                 const float& min_all_modified_ratio,
                                                                                 const float& min_overlap_modified_ratio) {
   vector<typename unordered_map<string, bitset<NUM_PROTEOFORMS>>::const_iterator> nav;
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

vector<string> convert(const unordered_set<string>& a_set) {
   vector<string> result;
   result.assign(a_set.begin(), a_set.end());
   sort(result.begin(), result.end());
   return result;
}

// entities_column argument starts counting at 0
index_to_entitites_phegen_result getIndexToEntititesPheGen(const string& path_file_PheGenI_full, const double& max_p_value) {
   ifstream file_phegen(path_file_PheGenI_full);
   string line, field, trait, gene, gene2;
   string p_value_str;
   long double p_value;
   unordered_set<string> temp_gene_set, temp_trait_set;
   vector<string> index_to_genes, index_to_traits;

   if (!file_phegen.is_open()) {
      throw runtime_error("Cannot open " + path_file_PheGenI_full);
   }

   getline(file_phegen, line);                  // Read header line
   while (getline(file_phegen, field, '\t')) {  // Read #
      getline(file_phegen, trait, '\t');        // Read Trait
      getline(file_phegen, field, '\t');        // Read SNP rs
      getline(file_phegen, field, '\t');        // Read Context
      getline(file_phegen, gene, '\t');         //	Gene
      getline(file_phegen, field, '\t');        //	Gene ID
      getline(file_phegen, gene2, '\t');        //	Gene 2
      getline(file_phegen, field, '\t');        //	Gene ID 2
      getline(file_phegen, field, '\t');        // Read Chromosome
      getline(file_phegen, field, '\t');        // Read Location
      getline(file_phegen, p_value_str, '\t');  // Read P-Value
      getline(file_phegen, line);               // Skip header line leftoever: Source,	PubMed,	Analysis ID,	Study ID,	Study Name

      try {
         p_value = stold(p_value_str);
         if (p_value <= GENOME_WIDE_SIGNIFICANCE) {
            temp_trait_set.insert(trait);
            temp_gene_set.insert(gene);
            temp_gene_set.insert(gene2);
         }
      } catch (const std::exception& ex) {
         cerr << "Error converting: **" << p_value_str << "**\n";
      }
   }

   //Create the final data structures
   index_to_genes = convert(temp_gene_set);
   index_to_traits = convert(temp_trait_set);

   return {index_to_genes, index_to_traits};
}

load_entitites_phegen_result loadEntitiesPheGen(const string& path_file_PheGenI_full, const double& max_p_value) {
   const auto [index_to_genes, index_to_traits] = getIndexToEntititesPheGen(path_file_PheGenI_full, max_p_value);
   const auto genes_to_index = getEntitiesToIndex(index_to_genes);
   const auto traits_to_index = getEntitiesToIndex(index_to_traits);

   cerr << "PHEGEN genes: " << index_to_genes.size() << " = " << genes_to_index.size() << "\n";
   cerr << "PHEGEN traits: " << index_to_traits.size() << " = " << traits_to_index.size() << "\n";

   return {index_to_genes, index_to_traits, genes_to_index, traits_to_index};
}

load_trait_gene_sets_result loadTraitGeneSets(const string& path_file_phegen,
                                              const double& max_p_value,
                                              const vector<string>& index_to_genes,
                                              const vector<string>& index_to_traits,
                                              const unordered_map<string, int>& genes_to_index,
                                              const unordered_map<string, int>& traits_to_index) {
   ifstream file_phegen(path_file_phegen);
   string line, field, trait, gene, gene2;
   string p_value_str;
   long double p_value;
   unordered_map<string, bitset<NUM_PHEGEN_GENES>> sets_to_genes;
   unordered_map<string, bitset<NUM_PHEGEN_TRAITS>> genes_to_sets;

   if (!file_phegen.is_open()) {
      throw runtime_error("Cannot open " + path_file_phegen);
   }

   getline(file_phegen, line);                  // Read header line
   while (getline(file_phegen, field, '\t')) {  // Read #
      getline(file_phegen, trait, '\t');        // Read Trait
      getline(file_phegen, field, '\t');        // Read SNP rs
      getline(file_phegen, field, '\t');        // Read Context
      getline(file_phegen, gene, '\t');         //	Gene
      getline(file_phegen, field, '\t');        //	Gene ID
      getline(file_phegen, gene2, '\t');        //	Gene 2
      getline(file_phegen, field, '\t');        //	Gene ID 2
      getline(file_phegen, field, '\t');        // Read Chromosome
      getline(file_phegen, field, '\t');        // Read Location
      getline(file_phegen, p_value_str, '\t');  // Read P-Value
      getline(file_phegen, line);               // Skip header line leftoever: Source,	PubMed,	Analysis ID,	Study ID,	Study Name

      try {
         p_value = stold(p_value_str);
         // cerr << "Converted correctly: " << p_value_str << "\n";
         if (p_value <= GENOME_WIDE_SIGNIFICANCE) {
            sets_to_genes[trait].set(genes_to_index.at(gene));
            sets_to_genes[trait].set(genes_to_index.at(gene2));
            genes_to_sets[gene].set(traits_to_index.at(trait));
            genes_to_sets[gene2].set(traits_to_index.at(trait));
         }
      } catch (const std::exception& ex) {
         cerr << "Error converting: **" << p_value_str << "**\n";
      }
   }

   return {genes_to_sets, sets_to_genes};
}

// Loads mapping from one column of strings to a second column of strings. If there are multiple values on the second column for one value on the left, then
//they are separated by spaces. Columns are separated by tabs.
load_mapping_result loadMapping(const string& path_file_mapping) {
   unordered_multimap<string, string> ones_to_others;
   unordered_multimap<string, string> others_to_ones;
   ifstream file_mapping(path_file_mapping);
   string to_str, from_str;

   if (!file_mapping.is_open()) {
      throw runtime_error(path_file_mapping);
   }

   getline(file_mapping, from_str);  // Discard the header line.
   while (getline(file_mapping, from_str, '\t')) {
      char c;

      while (file_mapping.get(c)) {
         if (c == ' ') {
            ones_to_others.emplace(from_str, to_str);
            others_to_ones.emplace(to_str, from_str);
            to_str = "";
         } else if (c == '\n') {
            ones_to_others.emplace(from_str, to_str);
            others_to_ones.emplace(to_str, from_str);
            to_str = "";
            break;
         } else {
            to_str = to_str.append(1, c);
         }
      }
   }
   cerr << "Mapping loaded.\n";
   return {ones_to_others, others_to_ones};
}

unordered_map<string, bitset<NUM_PHEGEN_PROTEINS>> convertGeneSetsToProteinSets(const unordered_map<string, bitset<NUM_PHEGEN_GENES>>& traits_to_genes,
                                                                                const vector<string>& index_to_genes,
                                                                                const unordered_multimap<string, string>& mapping_genes_to_proteins,
                                                                                const unordered_map<string, int>& proteins_to_index,
                                                                                const unordered_multimap<string, string>& adjacency_list_proteins) {
   unordered_map<string, bitset<NUM_PHEGEN_PROTEINS>> traits_to_proteins;
   for (const auto& trait_to_genes : traits_to_genes) {  // For each trait entry
      unordered_set<string> candidates;
      for (int I = 0; I < NUM_PHEGEN_GENES; I++) {  // For each gene in this trait
         if (trait_to_genes.second.test(I)) {
            auto range = mapping_genes_to_proteins.equal_range(index_to_genes[I]);  // For each mapping protein
            for (auto it = range.first; it != range.second; it++) {
               candidates.insert(it->second);
            }
         }
      }

      // Keep only those connected to any of the other gene set members in the reference protein network
      for (const auto& candidate : candidates) {
         auto range = adjacency_list_proteins.equal_range(candidate);
         for (auto it = range.first; it != range.second; ++it) {
            if (candidates.find(it->second) != candidates.end()) {
               traits_to_proteins[trait_to_genes.first].set(proteins_to_index.at(candidate));  // Set in stone that the candidate is in the new set
               break;
            }
         }
      }
   }

   return traits_to_proteins;
}

template <size_t total_num_result_entities, size_t total_num_prev_entities>
unordered_map<string, bitset<total_num_result_entities>> convertProteinSetsToProteoformSets(const map<string, bitset<total_num_prev_entities>>& traits_to_prev_entities,
                                                                                            const vector<string>& index_to_prev_entities,
                                                                                            const multimap<string, string>& mapping,
                                                                                            const map<string, int>& proteins_to_index,
                                                                                            const multimap<string, string>& adjacency_list_proteins) {
   unordered_map<string, bitset<total_num_result_entities>> traits_to_entities;
   for (const auto& trait_to_genes : traits_to_prev_entities) {  // For each trait entry
      unordered_set<string> candidates;
      for (int I = 0; I < NUM_PHEGEN_PROTEINS; I++) {  // For each member in this trait
         if (trait_to_genes.second.test(I)) {
            auto range = mapping.equal_range(index_to_prev_entities[I]);  // For each mapping entity
            for (auto it = range.first; it != range.second; it++) {
               candidates.insert(it->second);
            }
         }
      }

      // Keep only those connected to any of the other gene set members in the reference protein network
      for (const auto& candidate : candidates) {
         auto range = adjacency_list_proteins.equal_range(candidate);
         for (auto it = range.first; it != range.second; ++it) {
            if (candidates.find(it->second) != candidates.end()) {
               traits_to_entities[trait_to_genes.first].set(proteins_to_index.at(candidate));  // Set in stone that the candidate is in the new set
               break;
            }
         }
      }
   }

   return traits_to_entities;
}