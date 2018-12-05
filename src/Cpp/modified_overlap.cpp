#include "modified_overlap.hpp"

using namespace std;

namespace modified_overlap {

template <size_t total_proteoforms>
void writeReportRecords(ofstream& output,
                        const map<pair<string, string>, bitset<total_proteoforms>>& examples,
                        const unordered_map<string, string>& pathways_to_names,
                        const unordered_map<string, bitset<total_proteoforms>>& pathways_to_proteoforms,
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
                          const map<pair<string, string>, bitset<NUM_PHEGEN_PROTEOFORMS>>& examples,
                          const unordered_map<string, string>& pathways_to_names,
                          const unordered_map<string, bitset<NUM_PHEGEN_PROTEOFORMS>>& pathways_to_proteoforms,
                          const vector<string> index_to_proteoforms) {
   output << "PHENOTYPE_1\tPHENOTYPE_1_2\tPHENOTYPE_1_NAME\tPHENOTYPE_2_NAME\t";
   output << "PHENOTYPE_1_PROTEOFORM_SIZE\tPHENOTYPE_2_PROTEOFORM_SIZE\tOVERLAP_SIZE\t";
   output << "OVERLAP_PROTEOFORMS\n";
   writeReportRecords(output, examples, pathways_to_names, pathways_to_proteoforms, index_to_proteoforms);
}



void writeFrequencies(std::string_view modifications_file_path, std::string_view proteins_file_path, std::string_view proteoforms_file_path,
                      const map<string, int>& mod_to_freq, const map<string, int>& proteins_to_freq, const map<string, int>& proteoforms_to_freq) {
   ofstream modifications_file(modifications_file_path.data());
   ofstream proteins_file(proteins_file_path.data());
   ofstream proteoforms_file(proteoforms_file_path.data());

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

template <size_t total_proteoforms>
Frequencies getFrequencies(const map<pair<string, string>, bitset<total_proteoforms>>& proteoform_overlap_pairs,
                           const vector<string>& index_to_proteoforms) {
   Frequencies frequencies;
   for (const auto& overlap_pair : proteoform_overlap_pairs) {
      // For each overlap set
      for (int I = 0; I < overlap_pair.second.size(); I++) {
         if (overlap_pair.second.test(I)) {
            string accession = proteoform::getAccession(index_to_proteoforms[I]);
            if (frequencies.proteins.find(accession) == frequencies.proteins.end()) {
               frequencies.proteins.emplace(accession, 0);
            }
            frequencies.proteins[accession]++;

            if (frequencies.proteoforms.find(index_to_proteoforms[I]) == frequencies.proteoforms.end()) {
               frequencies.proteoforms.emplace(index_to_proteoforms[I], 0);
            }
            frequencies.proteoforms[index_to_proteoforms[I]]++;

            for (const auto& modification : proteoform::getModifications(index_to_proteoforms[I])) {
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

void plotFrequencies(const std::string& report_file_path, const std::string& modifications_file_path, const std::string& proteins_file_path, const std::string& proteoforms_file_path) {
   string command = "Rscript src/4_modified_overlap/modified_overlap.R " + report_file_path + " " +
                    modifications_file_path + " " + proteins_file_path + " " + proteoforms_file_path;
   cerr << "Plotting modified overlap frequencies: " << command << endl;
   system(command.c_str());
}

void reportPathwayPairs(std::string_view path_file_proteoform_search,
                        std::string_view report_file_path,
                        std::string_view modifications_file_path,
                        std::string_view proteins_file_path,
                        std::string_view proteoforms_file_path) {
   const auto [index_to_proteoforms, proteoforms_to_index] = pathway::readEntities(path_file_proteoform_search);
   const auto pathways_to_names = loadPathwayNames(path_file_proteoform_search);
   const auto pathways_to_proteoforms = loadProteoformSets(path_file_proteoform_search, proteoforms_to_index, true);
   const bitset<NUM_PROTEOFORMS> modified_proteoforms = proteoform::getSetOfModifiedProteoforms<NUM_PROTEOFORMS>(index_to_proteoforms);

   cout << "Reporting pathway pairs with only modified overlap...\n";

   // Compare all the pairs of selected pathways and select pairs that overlap only in a percentage of modified proteins
   const auto& examples = findOverlappingProteoformSets(pathways_to_proteoforms, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE,
                                                        modified_proteoforms, MIN_MODIFIED_ALL_MEMBERS_RATIO, MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

   ofstream report(report_file_path.data());
   writePathwayReport(report, examples, pathways_to_names, pathways_to_proteoforms, index_to_proteoforms);

   const auto [mod_to_freq, proteins_to_freq, proteoforms_to_freq] = getFrequencies(examples, index_to_proteoforms);
   writeFrequencies(modifications_file_path, proteins_file_path, proteoforms_file_path, mod_to_freq, proteins_to_freq, proteoforms_to_freq);
   plotFrequencies(report_file_path.data(), modifications_file_path.data(), proteins_file_path.data(), proteoforms_file_path.data());
}

void reportPhenotypePairs(std::string_view path_file_gene_search,
                          std::string_view path_file_protein_search,
                          std::string_view path_file_proteoform_search,
                          std::string_view path_file_PheGenI_full,
                          std::string_view path_file_mapping_proteins_genes,
                          std::string_view report_file_path,
                          std::string_view modifications_file_path,
                          std::string_view proteins_file_path,
                          std::string_view proteoforms_file_path) {
   // Load data: trait gene sets, trait proteoform sets
   cout << "Loading PheGen data\n";
   const auto [reactome_index_to_genes, reactome_genes_to_index] = pathway::readEntities(path_file_gene_search);
   const auto [index_to_genes, index_to_traits, genes_to_index, traits_to_index] = loadGenesPheGen(path_file_PheGenI_full, GENOME_WIDE_SIGNIFICANCE, reactome_genes_to_index);
   const auto [proteins_to_genes, genes_to_proteins] = loadMapping(path_file_mapping_proteins_genes.data());
   const auto [proteoforms_to_proteins, proteins_to_proteoforms] = loadMapping(path_file_proteoform_search.data());
   const auto [index_to_proteins, proteins_to_index] = deductProteinsFromGenes(path_file_mapping_proteins_genes, genes_to_index, genes_to_proteins);
   const auto [index_to_proteoforms, proteoforms_to_index] = deductProteoformsFromProteins(proteins_to_proteoforms, proteins_to_index);
   const bitset<NUM_PHEGEN_PROTEOFORMS> modified_proteoforms = proteoform::getSetOfModifiedProteoforms<NUM_PHEGEN_PROTEOFORMS>(index_to_proteoforms);

   const auto [adjacency_list_proteins, adjacency_list_proteoforms] = loadReactomeNetworks(path_file_gene_search, path_file_protein_search, path_file_proteoform_search);
   const auto [genes_to_traits, traits_to_genes] = loadTraitGeneSets(path_file_PheGenI_full.data(), GENOME_WIDE_SIGNIFICANCE, index_to_genes, index_to_traits, genes_to_index, traits_to_index, reactome_genes_to_index);
   const auto sets_to_names = createTraitNames(traits_to_genes);
   const unordered_map<string, bitset<NUM_PHEGEN_PROTEINS>> traits_to_proteins = convertGeneSets(traits_to_genes, index_to_genes, genes_to_proteins, proteins_to_index, adjacency_list_proteins);
   const unordered_map<string, bitset<NUM_PHEGEN_PROTEOFORMS>> traits_to_proteoforms = convertProteinSets(traits_to_proteins, index_to_proteins, proteins_to_proteoforms, proteoforms_to_index, adjacency_list_proteoforms);

   // Calculate overlap
   const auto overlapping_proteoform_set_pairs = findOverlappingPairs(traits_to_proteoforms, 1, NUM_PHEGEN_PROTEOFORMS, 1, NUM_PHEGEN_PROTEOFORMS);
   const auto non_overlapping_proteoform_set_pairs = findOverlappingPairs(traits_to_proteoforms, 0, 0, 1, NUM_PHEGEN_PROTEOFORMS);
   const auto& examples = findOverlappingProteoformSets(traits_to_proteoforms, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE,
                                                        modified_proteoforms, MIN_MODIFIED_ALL_MEMBERS_RATIO, MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

   // Write report
   ofstream report(report_file_path.data());
   writePhenotypeReport(report, examples, sets_to_names, traits_to_proteoforms, index_to_proteoforms);

   // Analyse modifications
   const auto [mod_to_freq, proteins_to_freq, proteoforms_to_freq] = getFrequencies(examples, index_to_proteoforms);
   writeFrequencies(modifications_file_path, proteins_file_path, proteoforms_file_path, mod_to_freq, proteins_to_freq, proteoforms_to_freq);
   plotFrequencies(report_file_path.data(), modifications_file_path.data(), proteins_file_path.data(), proteoforms_file_path.data());
}

void doAnalysis(std::string_view path_file_gene_search,
                std::string_view path_file_protein_search,
                std::string_view path_file_proteoform_search,
                std::string_view path_file_PheGenI_full,
                std::string_view path_file_mapping_proteins_to_genes,
                std::string_view path_file_report_pathway,
                std::string_view path_file_modified_overlap_pathway_proteins,
                std::string_view path_file_modified_overlap_pathway_proteoforms,
                std::string_view path_file_modified_overlap_pathway_modifications,
                std::string_view path_file_report_trait,
                std::string_view path_file_modified_overlap_trait_modifications,
                std::string_view path_file_modified_overlap_trait_proteins,
                std::string_view path_file_modified_overlap_trait_proteoforms) {
   cout << "Searching for modified overlap examples...\n";

   // Part 1: Find examples of pathways that overlap only in modified proteins
   reportPathwayPairs(path_file_proteoform_search,
                      path_file_report_pathway,
                      path_file_modified_overlap_pathway_modifications,
                      path_file_modified_overlap_pathway_proteins,
                      path_file_modified_overlap_pathway_proteoforms);

   // Part 2: Find examples of disease modules that overlap only in modified proteins
   reportPhenotypePairs(path_file_gene_search,
                        path_file_protein_search,
                        path_file_proteoform_search,
                        path_file_PheGenI_full,
                        path_file_mapping_proteins_to_genes,
                        path_file_report_trait,
                        path_file_modified_overlap_trait_modifications,
                        path_file_modified_overlap_trait_proteins,
                        path_file_modified_overlap_trait_proteoforms);
}

}  // namespace modified_overlap
