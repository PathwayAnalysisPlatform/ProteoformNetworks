#include "dataset.hpp"

namespace pathway {

dataset::dataset(std::string_view path_file_gene_mapping,
                 std::string_view path_file_protein_mapping,
                 std::string_view path_file_proteoform_mapping) {
   // Create pathway sets
   setPathwayNames(path_file_protein_mapping);
   setGeneMapping(path_file_gene_mapping);
   setProteinMapping(path_file_protein_mapping);
   setProteoformMapping(path_file_proteoform_mapping);

   // Create interaction network
   // calculateInteractionNetworks();
}

std::string dataset::getName() const {
   return name;
}

void dataset::setName(std::string_view value) {
   name = value;
}

const std::unordered_map<std::string, std::string>& dataset::getPathwayNames() const {
   return pathways_to_names;
}

const entities_bimap& dataset::getGenes() const {
   return genes;
}
const entities_bimap& dataset::getProteins() const {
   return proteins;
}
const entities_bimap& dataset::getProteoforms() const {
   return proteoforms;
}

const std::unordered_map<std::string, std::unordered_set<std::string>>& dataset::getGenesToReactions() const {
   return genes_to_reactions;
}
const std::unordered_map<std::string, std::unordered_set<std::string>>& dataset::getGenesToPathways() const {
   return genes_to_pathways;
}
const std::unordered_map<std::string, std::unordered_set<std::string>>& dataset::getProteinsToReactions() const {
   return proteins_to_reactions;
}
const std::unordered_map<std::string, std::unordered_set<std::string>>& dataset::getProteinsToPathways() const {
   return proteins_to_pathways;
}
const std::unordered_map<std::string, std::unordered_set<std::string>>& dataset::getProteoformsToReactions() const {
   return proteoforms_to_reactions;
}
const std::unordered_map<std::string, std::unordered_set<std::string>>& dataset::getProteoformsToPathways() const {
   return proteoforms_to_pathways;
}

void dataset::setPathwayNames(std::string_view path_file_mapping) {
   std::cerr << "Setting pathway names\n";
   std::ifstream file_search(path_file_mapping.data());
   std::string field, pathway_id, pathway_name, lefover;

   getline(file_search, lefover);              // Skip csv header line
   while (getline(file_search, field, ',')) {  // Read PROTEIN
      getline(file_search, field, ',');        // Read REACTION_STID
      getline(file_search, pathway_id, ',');   // Read PATHWAY_STID
      getline(file_search, field, ',');        // Read REACTION_NAME
      getline(file_search, pathway_name);      // Read PATHWAY_NAME
      pathways_to_names.emplace(pathway_id, pathway_name);
   }

   std::cerr << "Finished setting pathway names\n";
}

void dataset::setGeneMapping(std::string_view path_file_mapping) {
   genes = readEntities(path_file_mapping);
   std::cerr << "Read " << genes.index_to_entities.size() << " genes.\n";

   std::ifstream map_file(path_file_mapping.data());
   std::string protein, gene, reaction, pathway, leftover;

   if (!map_file.is_open()) {
      throw std::runtime_error("Could not open gene mapping file.");
   }

   getline(map_file, leftover);            // Skip csv header line
   while (getline(map_file, gene, ',')) {  // Read GENE
      getline(map_file, reaction, ',');    // Read REACTION_STID
      getline(map_file, pathway, ',');     // Read PATHWAY_STID
      getline(map_file, protein);          // Read PROTEIN

      pathways_to_genes[pathway].set(genes.entities_to_index.at(gene));
      genes_to_pathways[gene].insert(pathway);

      reactions_to_genes[reaction].set(genes.entities_to_index.at(gene));
      genes_to_reactions[gene].insert(reaction);
   }
   std::cerr << "Filled " << pathways_to_genes.size() << " pathways for genes.\n";
   std::cerr << "Filled " << reactions_to_genes.size() << " reactions for genes.\n";
}

void dataset::setProteinMapping(std::string_view path_file_mapping) {
   proteins = readEntities(path_file_mapping.data());
   std::cerr << "Read " << proteins.index_to_entities.size() << " proteins.\n";

   std::ifstream map_file(path_file_mapping.data());
   std::string protein, reaction, pathway, leftover;

   if (!map_file.is_open())
      throw std::runtime_error("Could not open protein mapping file.");

   getline(map_file, leftover);               // Skip csv header line
   while (getline(map_file, protein, ',')) {  // Read PROTEIN
      getline(map_file, reaction, ',');       // Read REACTION_STID
      getline(map_file, pathway, ',');        // Read PATHWAY_STID
      getline(map_file, leftover);            // Read rest of line

      pathways_to_proteins[pathway].set(proteins.entities_to_index.at(protein));
      proteins_to_pathways[protein].insert(pathway);

      reactions_to_proteins[reaction].set(proteins.entities_to_index.at(protein));
      proteins_to_reactions[protein].insert(reaction);
   }
   std::cerr << "Filled " << pathways_to_proteins.size() << " pathways for proteins.\n";
   std::cerr << "Filled " << reactions_to_proteins.size() << " reactions for proteins.\n";
}

void dataset::setProteoformMapping(std::string_view path_file_mapping) {
   proteoforms = readEntities(path_file_mapping);
   std::cerr << "Read " << proteoforms.index_to_entities.size() << " proteoforms.\n";

   std::ifstream map_file(path_file_mapping.data());
   std::string proteoform, reaction, pathway, leftover;

   if (!map_file.is_open())
      throw std::runtime_error("Could not open proteoform mapping file.");

   getline(map_file, leftover);  // Skip csv header line
   while (map_file.peek() != EOF) {
      if (map_file.peek() == '\"')  // Read initial "
         map_file.get();
      map_file.get();  // Read initial " or [
      getline(map_file, proteoform, ']');
      getline(map_file, leftover, ',');  // Read end " if it exist
      getline(map_file, reaction, ',');  // Read REACTION_STID
      getline(map_file, pathway);        // Read PATHWAY_STID

      pathways_to_proteoforms[pathway].set(proteoforms.entities_to_index.at(proteoform));
      proteoforms_to_pathways[proteoform].insert(pathway);

      reactions_to_proteoforms[reaction].set(proteoforms.entities_to_index.at(proteoform));
      proteoforms_to_reactions[proteoform].insert(reaction);
   }
   std::cerr << "Filled " << pathways_to_proteoforms.size() << " pathways for proteoforms.\n";
   std::cerr << "Filled " << reactions_to_proteoforms.size() << " reactions for proteoforms.\n";
}

void dataset::calculateInteractionNetworks() {
   throw std::runtime_error("Not implemented exception");
}

}  // namespace pathway
