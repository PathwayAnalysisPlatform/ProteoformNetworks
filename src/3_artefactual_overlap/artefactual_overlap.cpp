#include "artefactual_overlap.hpp"

using namespace std;

namespace artefactual_overlap {

/**
 * Find artefactual overlaps: pairs of pathways that share nodes only in the
 * gene or protein level, but not at the proteoform level
 */
template <size_t size_genes, size_t size_proteins, size_t size_proteoforms>
set<pair<string, string>> findArtefactualOverlapPairs(const unordered_map<string, bitset<size_proteoforms>>& sets_to_proteoforms,
                                                      const map<pair<string, string>, bitset<size_genes>>& overlapping_gene_set_pairs,
                                                      const map<pair<string, string>, bitset<size_proteins>>& overlapping_protein_set_pairs,
                                                      const map<pair<string, string>, bitset<size_proteoforms>>& non_overlapping_proteoform_set_pairs) {
   set<pair<string, string>> result;

   cerr << "Comparing gene and proteoform pairs..." << endl;
   for (const auto& gene_set_pair : overlapping_gene_set_pairs) {                                                          // For each pair of overlapping sets at gene level
      if (non_overlapping_proteoform_set_pairs.find(gene_set_pair.first) != non_overlapping_proteoform_set_pairs.end()) {  // Check if they do not overlap
         result.insert(gene_set_pair.first);
      }
   }
   cerr << "Comparing protein and proteoform pairs..." << endl;
   for (const auto& protein_set_pair : overlapping_protein_set_pairs) {
      if (non_overlapping_proteoform_set_pairs.find(protein_set_pair.first) != non_overlapping_proteoform_set_pairs.end()) {
         result.insert(protein_set_pair.first);
      }
   }
   cout << "Found " << result.size() << " examples.\n";
   return result;
}

template <size_t size_proteoforms>
bitset<size_proteoforms> getProteoformsWithAccessions(const unordered_set<string>& accessions,
                                                      const bitset<size_proteoforms>& proteoform_set,
                                                      const vector<string>& index_to_proteoforms) {
   bitset<size_proteoforms> result;
   for (int I = 0; I < proteoform_set.size(); I++) {
      if (proteoform_set.test(I)) {
         if (accessions.find(getAccession(index_to_proteoforms.at(I))) != accessions.end()) {
            result.set(I);
         }
      }
   }
   return result;
}

template <size_t size_genes, size_t size_proteins, size_t size_proteoforms>
void writeReportRecords(ofstream& output,
                        const set<pair<string, string>> examples,
                        const unordered_map<string, string>& sets_to_names,
                        const unordered_map<string, bitset<size_genes>>& sets_to_genes,
                        const unordered_map<string, bitset<size_proteins>>& sets_to_proteins,
                        const unordered_map<string, bitset<size_proteoforms>>& sets_to_proteoforms,
                        const vector<string> index_to_genes,
                        const vector<string> index_to_proteins,
                        const vector<string> index_to_proteoforms) {
   cout << "Writing records...\n";

   auto t0 = clock();

   for (const auto& example : examples) {
      bitset<size_genes> overlap_genes = sets_to_genes.at(example.first) & sets_to_genes.at(example.second);
      bitset<size_proteins> overlap_proteins = sets_to_proteins.at(example.first) & sets_to_proteins.at(example.second);
      bitset<size_proteoforms> overlap_proteoforms = sets_to_proteoforms.at(example.first) & sets_to_proteoforms.at(example.second);
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
      printMembers(output, overlap_genes, index_to_genes);
      output << "\t";
      printMembers(output, overlap_proteins, index_to_proteins);
      output << "\t";
      printMembers(output, overlap_proteoforms, index_to_proteoforms);
      output << "\t";

      printMembers(output, decomposed_overlap_proteoforms_1, index_to_proteoforms);
      output << "\t";
      printMembers(output, decomposed_overlap_proteoforms_2, index_to_proteoforms);
      output << "\n";
   }
   auto t1 = clock();
   cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
}

void writePathwayReport(ofstream& output,
                        const set<pair<string, string>> examples,
                        const unordered_map<string, string>& sets_to_names,
                        const unordered_map<string, bitset<NUM_GENES>>& sets_to_genes,
                        const unordered_map<string, bitset<NUM_PROTEINS>>& sets_to_proteins,
                        const unordered_map<string, bitset<NUM_PROTEOFORMS>>& sets_to_proteoforms,
                        const vector<string> index_to_genes,
                        const vector<string> index_to_proteins,
                        const vector<string> index_to_proteoforms) {
   cerr << "Reporting pathway pairs with artefactual overlap...\n";
   output << "PATHWAY_1\tPATHWAY_2\tPATHWAY_1_NAME\tPATHWAY_2_NAME\t";
   output << "PATHWAY_1_GENE_SIZE\tPATHWAY_1_PROTEIN_SIZE\tPATHWAY_1_PROTEOFORM_SIZE\t";
   output << "PATHWAY_2_GENE_SIZE\tPATHWAY_2_PROTEIN_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\t";
   output << "GENE_OVERLAP\tPROTEIN_OVERLAP\tPROTEOFORM_OVERLAP\t";
   output << "OVERLAP_GENES\tOVERLAP_PROTEINS\tOVERLAP_PROTEOFORMS\tDECOMPOSED_OVERLAP_PROTEOFORMS_1\tDECOMPOSED_OVERLAP_PROTEOFORMS_2\n";
   writeReportRecords(output, examples, sets_to_names,
                      sets_to_genes, sets_to_proteins, sets_to_proteoforms,
                      index_to_genes, index_to_proteins, index_to_proteoforms);
}

void writePhenotypeReport(ofstream& output,
                          const set<pair<string, string>> examples,
                          const unordered_map<string, string>& sets_to_names,
                          const unordered_map<string, bitset<NUM_PHEGEN_GENES>>& sets_to_genes,
                          const unordered_map<string, bitset<NUM_PHEGEN_PROTEINS>>& sets_to_proteins,
                          const unordered_map<string, bitset<NUM_PHEGEN_PROTEOFORMS>>& sets_to_proteoforms,
                          const vector<string> index_to_genes,
                          const vector<string> index_to_proteins,
                          const vector<string> index_to_proteoforms) {
   output << "PHENOTYPE_1\tPHENOTYPE_2\tPHENOTYPE_1_NAME\tPHENOTYPE_2_NAME\t";
   output << "PHENOTYPE_1_GENE_SIZE\tPHENOTYPE_1_PROTEIN_SIZE\tPHENOTYPE_1_PROTEOFORM_SIZE\t";
   output << "PHENOTYPE_2_GENE_SIZE\tPHENOTYPE_2_PROTEIN_SIZE\tPHENOTYPE_2_PROTEOFORM_SIZE\t";
   output << "GENE_OVERLAP\tPROTEIN_OVERLAP\tPROTEOFORM_OVERLAP\t";
   output << "OVERLAP_GENES\tOVERLAP_PROTEINS\tOVERLAP_PROTEOFORMS\tDECOMPOSED_OVERLAP_PROTEOFORMS_1\tDECOMPOSED_OVERLAP_PROTEOFORMS_2\n";

   writeReportRecords(output, examples, sets_to_names,
                      sets_to_genes, sets_to_proteins, sets_to_proteoforms,
                      index_to_genes, index_to_proteins, index_to_proteoforms);
}

template <size_t total_num_entities>
void writeSetReport(const string& path_file_report,
                    const map<string, bitset<total_num_entities>>& sets_to_entities,
                    const vector<string>& index_to_entities) {
   ofstream file_report(path_file_report);
   cout << "Writing report: " << path_file_report << "\n";
   if (!file_report.is_open()) {
      throw runtime_error("Cannot write to " + path_file_report);
   }

   for (const auto& set_to_entities : sets_to_entities) {
      file_report << set_to_entities.first << "\t";
      for (int I = 0; I < total_num_entities; I++) {
         if (set_to_entities.second.test(I)) {
            file_report << index_to_entities.at(I) << "\t";
         }
      }
      file_report << "\n";
   }
}

void writeEntitiesReport(const string& path_file_report, const vector<string>& entities) {
   ofstream file_report(path_file_report);

   cerr << "Writing report: " << path_file_report << "\n";
   if (!file_report.is_open()) {
      throw runtime_error("Cannot write to " + path_file_report);
   }

   for (const auto& entity : entities) {
      file_report << entity << "\n";
   }
}

Entities_bimap deductProteinsFromGenes(const string& path_file_mapping_proteins_genes,
                                       const unordered_map<string, int>& genes_to_index,
                                       const unordered_multimap<string, string>& genes_to_proteins) {
   unordered_set<string> temp_protein_set;

   for (const auto& gene_entry : genes_to_index) {
      auto ret = genes_to_proteins.equal_range(gene_entry.first);
      for (auto it = ret.first; it != ret.second; ++it) {
         temp_protein_set.insert(it->second);
      }
   }

   vector<string> index_to_proteins = convert(temp_protein_set);
   unordered_map<string, int> proteins_to_index = getEntitiesToIndex(index_to_proteins);

   cerr << "PHEGEN proteins: " << index_to_proteins.size() << " = " << proteins_to_index.size() << "\n";

   return {index_to_proteins, proteins_to_index};
}

Entities_bimap deductProteoformsFromProteins(const unordered_multimap<string, string>& proteins_to_proteoforms, const unordered_map<string, int>& proteins_to_index) {
   unordered_set<string> temp_proteoform_set;

   for (const auto& entry : proteins_to_index) {
      auto ret = proteins_to_proteoforms.equal_range(entry.first);
      if (ret.first != ret.second) {
         for (auto it = ret.first; it != ret.second; ++it) {
            temp_proteoform_set.insert(it->second);
         }
      } else {
         temp_proteoform_set.insert(entry.first + ";");
      }
   }

   vector<string> index_to_proteoforms = convert(temp_proteoform_set);
   unordered_map<string, int> proteoforms_to_index = getEntitiesToIndex(index_to_proteoforms);

   cerr << "PHEGEN proteoforms: " << index_to_proteoforms.size() << " = " << proteoforms_to_index.size() << "\n";

   return {index_to_proteoforms, proteoforms_to_index};
}

map<string, bitset<NUM_GENES>> loadReactionsGeneMembers(const string& file_path, const map<string, int>& entities_to_index) {
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

unordered_multimap<string, string> loadGenesAdjacencyList(const string& search_file_path) {
   unordered_multimap<string, string> adjacenty_list;
   const auto [index_to_entities, entities_to_index] = loadEntities(search_file_path);
   const auto reactions_to_entities = loadGeneSets(search_file_path, entities_to_index, false);

   for (const auto& reaction_entry : reactions_to_entities) {
      vector<string> members;
      for (int I = 0; I < reaction_entry.second.size(); I++) {
         if (reaction_entry.second.test(I)) {
            members.push_back(index_to_entities[I]);
         }
      }

      for (const auto& one_member : members) {
         for (const auto& other_member : members) {
            adjacenty_list.emplace(one_member, other_member);
         }
      }
   }

   return adjacenty_list;
}

unordered_multimap<string, string> loadProteinsAdjacencyList(const string& search_file_path) {
   unordered_multimap<string, string> adjacenty_list;
   const auto [index_to_entities, entities_to_index] = loadEntities(search_file_path);
   const auto reactions_to_entities = loadProteinSets(search_file_path, entities_to_index, false);

   for (const auto& reaction_entry : reactions_to_entities) {
      vector<string> members;
      for (int I = 0; I < reaction_entry.second.size(); I++) {
         if (reaction_entry.second.test(I)) {
            members.push_back(index_to_entities[I]);
         }
      }

      for (const auto& one_member : members) {
         for (const auto& other_member : members) {
            adjacenty_list.emplace(one_member, other_member);
         }
      }
   }

   return adjacenty_list;
}

unordered_multimap<string, string> loadProteoformsAdjacencyList(const string& search_file_path) {
   unordered_multimap<string, string> adjacenty_list;
   const auto [index_to_entities, entities_to_index] = loadEntities(search_file_path);
   const auto reactions_to_entities = loadProteoformSets(search_file_path, entities_to_index, false);

   for (const auto& reaction_entry : reactions_to_entities) {
      vector<string> members;
      for (int I = 0; I < reaction_entry.second.size(); I++) {
         if (reaction_entry.second.test(I)) {
            members.push_back(index_to_entities[I]);
         }
      }

      for (const auto& one_member : members) {
         for (const auto& other_member : members) {
            adjacenty_list.emplace(one_member, other_member);
         }
      }
   }

   return adjacenty_list;
}

void reportPathwayPairs(const string& path_file_gene_search,
                        const string& path_file_protein_search,
                        const string& path_file_proteoform_search,
                        const string& report_file_path) {
   const auto [index_to_genes, genes_to_index] = loadEntities(path_file_gene_search);
   const auto [index_to_proteins, proteins_to_index] = loadEntities(path_file_protein_search);
   const auto [index_to_proteoforms, proteoforms_to_index] = loadEntities(path_file_proteoform_search);
   const auto pathways_to_names = loadPathwayNames(path_file_proteoform_search);

   unordered_map<string, bitset<NUM_GENES>> pathways_to_genes = loadGeneSets(path_file_gene_search, genes_to_index, true);
   unordered_map<string, bitset<NUM_PROTEINS>> pathways_to_proteins = loadProteinSets(path_file_protein_search, proteins_to_index, true);
   unordered_map<string, bitset<NUM_PROTEOFORMS>> pathways_to_proteoforms = loadProteoformSets(path_file_proteoform_search, proteoforms_to_index, true);

   cout << "Calculating gene sets overlap..." << endl;
   const auto overlapping_gene_set_pairs = findOverlappingPairs(pathways_to_genes, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE);
   cout << "Calculating protein sets overlap..." << endl;
   const auto overlapping_protein_set_pairs = findOverlappingPairs(pathways_to_proteins, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE);
   cout << "Calculating proteoform sets overlap..." << endl;
   const auto overlapping_proteoform_set_pairs = findOverlappingPairs(pathways_to_proteoforms, 1, NUM_PROTEOFORMS, 1, NUM_PROTEOFORMS);
   const auto non_overlapping_proteoform_set_pairs = findOverlappingPairs(pathways_to_proteoforms, 0, 0, 1, NUM_PROTEOFORMS);

   cerr << "Total proteoforms sets: " << pathways_to_proteoforms.size() << "\n";
   cerr << "Parejas: " << (pathways_to_proteoforms.size() * pathways_to_proteoforms.size() - pathways_to_proteoforms.size()) / 2 << "\n";
   cerr << "Overlapping: " << overlapping_proteoform_set_pairs.size() << "\n";
   cerr << "Non-overlapping: " << non_overlapping_proteoform_set_pairs.size() << "\n";
   cerr << "Total pairs defined: " << overlapping_proteoform_set_pairs.size() + non_overlapping_proteoform_set_pairs.size() << "\n";

   cout << "Finding examples of artifactual overlap...\n";
   const auto examples = findArtefactualOverlapPairs(pathways_to_proteoforms,
                                                     overlapping_gene_set_pairs,
                                                     overlapping_protein_set_pairs,
                                                     non_overlapping_proteoform_set_pairs);

   ofstream report(report_file_path);
   writePathwayReport(report, examples, pathways_to_names,
                      pathways_to_genes, pathways_to_proteins, pathways_to_proteoforms,
                      index_to_genes, index_to_proteins, index_to_proteoforms);
}

void reportPhenotypePairs(const std::string& path_file_gene_search,
                          const std::string& path_file_protein_search,
                          const std::string& path_file_proteoform_search,
                          const std::string& path_file_PheGenI_full,
                          const std::string& path_file_mapping_proteins_genes,
                          const std::string& path_file_report_trait) {
   // Load reference network

   // Load phegen sets

   // Convert to protein sets: reference gene network, phegen gene sets, mapping from genes to proteins
   // For each protein, decide if it is connected to the other proteins in the set or not in the reference network
   // == For each gene in the phegen set:
   // == Use the mapping from genes to proteins and find the corresponding proteins.
   // == For each correspoding protein:
   // == Add to the new phegen protein set if at least one of its neighbours in the reference network is in the phegen set
   // Convert to proteoform sets: reference protein network, phegen protein sets, mapping proteins to proteoforms
   // Find overlapping pairs
   // Write report

   // Data structures:
   // reference network: adjacency list: map<string, string>
   // phegen set: bitset<MAX_ENTITIES> members, index_to_entities, string set_name

   cout << "Loading PheGen data\n";
   const auto [reactome_index_to_genes, reactome_genes_to_index] = loadEntities(path_file_gene_search);
   const auto [index_to_genes, index_to_traits, genes_to_index, traits_to_index] = loadGenesPheGen(path_file_PheGenI_full, GENOME_WIDE_SIGNIFICANCE, reactome_genes_to_index);
   const auto [proteins_to_genes, genes_to_proteins] = loadMapping(path_file_mapping_proteins_genes);
   const auto [proteoforms_to_proteins, proteins_to_proteoforms] = loadMapping(path_file_proteoform_search);
   const auto [index_to_proteins, proteins_to_index] = deductProteinsFromGenes(path_file_mapping_proteins_genes, genes_to_index, genes_to_proteins);
   const auto [index_to_proteoforms, proteoforms_to_index] = deductProteoformsFromProteins(proteins_to_proteoforms, proteins_to_index);

   // Calculate adjacency lists for genes, proteins and proteoforms according to Reactome
   // An entity is neighbour of another gene if the participate in the same Reaction
   cout << "Loading Reactome network.\n";
   const unordered_multimap<string, string> adjacency_list_genes = loadGenesAdjacencyList(path_file_gene_search);
   const unordered_multimap<string, string> adjacency_list_proteins = loadProteinsAdjacencyList(path_file_protein_search);
   const unordered_multimap<string, string> adjacency_list_proteoforms = loadProteoformsAdjacencyList(path_file_proteoform_search);
   const auto [genes_to_traits, traits_to_genes] = loadTraitGeneSets(path_file_PheGenI_full, GENOME_WIDE_SIGNIFICANCE, index_to_genes, index_to_traits, genes_to_index, traits_to_index, reactome_genes_to_index);
   cerr << "Loaded gene sets\n";
   const unordered_map<string, bitset<NUM_PHEGEN_PROTEINS>> traits_to_proteins = convertGeneSets(traits_to_genes, index_to_genes, genes_to_proteins, proteins_to_index, adjacency_list_proteins);
   printMembers(cout, traits_to_genes.at("Abdominal Fat"), index_to_genes);
   cout << endl;
   printMembers(cout, traits_to_proteins.at("Abdominal Fat"), index_to_proteins);
   cout << endl;
   const unordered_map<string, bitset<NUM_PHEGEN_PROTEOFORMS>> traits_to_proteoforms = convertProteinSets(traits_to_proteins, index_to_proteins, proteins_to_proteoforms, proteoforms_to_index, adjacency_list_proteoforms);
   cout << "Finished calculating sets.\n";

   const auto overlapping_gene_set_pairs = findOverlappingPairs(traits_to_genes, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE);
   cout << "Calculating protein sets overlap..." << endl;
   const auto overlapping_protein_set_pairs = findOverlappingPairs(traits_to_proteins, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE);
   cout << "Calculating proteoform sets overlap..." << endl;
   const auto overlapping_proteoform_set_pairs = findOverlappingPairs(traits_to_proteoforms, 1, NUM_PHEGEN_PROTEOFORMS, 1, NUM_PHEGEN_PROTEOFORMS);
   const auto non_overlapping_proteoform_set_pairs = findOverlappingPairs(traits_to_proteoforms, 0, 0, 1, NUM_PHEGEN_PROTEOFORMS);

   cerr << "Total proteoforms sets: " << traits_to_proteoforms.size() << "\n";
   cerr << "Parejas: " << (traits_to_proteoforms.size() * traits_to_proteoforms.size() - traits_to_proteoforms.size()) / 2 << "\n";
   cerr << "Overlapping: " << overlapping_proteoform_set_pairs.size() << "\n";
   cerr << "Non-overlapping: " << non_overlapping_proteoform_set_pairs.size() << "\n";
   cerr << "Total pairs defined: " << overlapping_proteoform_set_pairs.size() + non_overlapping_proteoform_set_pairs.size() << "\n";

   cout << "Finding examples of artifactual overlap...\n";
   const auto examples = findArtefactualOverlapPairs(traits_to_proteoforms,
                                                     overlapping_gene_set_pairs,
                                                     overlapping_protein_set_pairs,
                                                     non_overlapping_proteoform_set_pairs);

   cout << "Writing report...\n";
   ofstream report(path_file_report_trait);

   unordered_map<string, string> sets_to_names;
   for (const auto& trait_entry : traits_to_genes) {
      sets_to_names.emplace(trait_entry.first, trait_entry.first);
   }
   writePhenotypeReport(report, examples, sets_to_names, traits_to_genes, traits_to_proteins, traits_to_proteoforms,
                        index_to_genes, index_to_proteins, index_to_proteoforms);
}

// Find pairs of modules/pathways that overlap on gene or protein network, but not in proteoform network
void doAnalysis(const std::string& path_file_gene_search,
                const std::string& path_file_protein_search,
                const std::string& path_file_proteoform_search,
                const std::string& path_file_PheGenI_full,
                const std::string& path_file_mapping_proteins_to_genes,
                const std::string& path_file_report_pathway,
                const std::string& path_file_report_trait) {
   cout << "Searching for artefactual overlap examples..." << endl;

   reportPathwayPairs(path_file_gene_search,
                      path_file_protein_search,
                      path_file_proteoform_search,
                      path_file_report_pathway);
   cout << "\n\n************************************\n\n";
   reportPhenotypePairs(path_file_gene_search,
                        path_file_protein_search,
                        path_file_proteoform_search,
                        path_file_PheGenI_full,
                        path_file_mapping_proteins_to_genes,
                        path_file_report_trait);
}

}  // namespace artefactual_overlap
