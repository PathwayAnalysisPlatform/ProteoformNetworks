#include "dataset.hpp"

namespace pathway {

dataset::dataset(std::string_view path_file_gene_mapping,
                 std::string_view path_file_protein_mapping,
                 std::string_view path_file_proteoform_mapping) {
   // Create pathway sets
   setGeneSets(path_file_gene_mapping);
   setProteinSets(path_file_protein_mapping);
   setProteoformSets(path_file_proteoform_mapping);
   setPathwayNames(path_file_protein_mapping);

   // Create interaction network
   calculateInteractionNetworks();
}

void dataset::setGeneSets(std::string_view path_file_mapping) {
   genes = readEntities(path_file_mapping);

   std::ifstream file_search(path_file_mapping.data());
   std::string field, gene, reaction, pathway;

   getline(file_search, field);               // Skip csv header line
   while (getline(file_search, gene, ',')) {  // Read GENE
      getline(file_search, field, ',');       // Read PROTEIN
      getline(file_search, reaction, ',');    // Read REACTION_STID
      getline(file_search, pathway, ',');     // Read PATHWAY_STID

      // result[pathways ? pathway : reaction].set(entities_to_index.find(gene)->second);
   }
}

void dataset::setPathwayNames(std::string_view path_file_mapping) {
}

void dataset::setProteinSets(std::string_view path_file_mapping) {
   proteins = readEntities(path_file_mapping.data());
}

void dataset::setProteoformSets(std::string_view path_file_mapping) {
   proteoforms = readEntities(path_file_mapping);
}

void dataset::calculateInteractionNetworks() {
   throw std::runtime_error("Not implemented exception");
}

}  // namespace pathway