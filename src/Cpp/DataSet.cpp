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
   std::cerr << "Calculating networks...\n";
   calculateInteractionNetworks();
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

const int dataset::getNumGenes() const {
   return genes.index_to_entities.size();
}
const int dataset::getNumProteins() const {
   return proteins.index_to_entities.size();
}
const int dataset::getNumProteoforms() const {
   return proteoforms.index_to_entities.size();
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
const std::unordered_multimap<std::string, std::string>& dataset::getGenesToProteins() const {
   return genes_to_proteins;
}
const std::unordered_multimap<std::string, std::string>& dataset::getProteinsToProteoforms() const {
   return proteins_to_proteoforms;
}
const std::unordered_multimap<std::string, std::string>& dataset::getGeneNetwork() const {
   return gene_network;
}
const std::unordered_multimap<std::string, std::string>& dataset::getProteinNetwork() const {
   return protein_network;
}
const std::unordered_multimap<std::string, std::string>& dataset::getProteoformNetwork() const {
   return proteoform_network;
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
   std::cerr << "Read " << getNumGenes() << " genes.\n";

   std::ifstream map_file(path_file_mapping.data());
   std::string protein, gene, reaction, pathway, leftover;
   std::set<std::pair<std::string, std::string>> temp_genes_to_proteins;

   if (!map_file.is_open()) {
      throw std::runtime_error("Could not open gene mapping file.");
   }

   getline(map_file, leftover);            // Skip csv header line
   while (getline(map_file, gene, ',')) {  // Read GENE
      getline(map_file, reaction, ',');    // Read REACTION_STID
      getline(map_file, pathway, ',');     // Read PATHWAY_STID
      getline(map_file, protein);          // Read PROTEIN

      temp_genes_to_proteins.insert(std::make_pair(gene, protein));

      pathways_to_genes[pathway].set(genes.entities_to_index.at(gene));
      genes_to_pathways[gene].insert(pathway);

      reactions_to_genes[reaction].set(genes.entities_to_index.at(gene));
      genes_to_reactions[gene].insert(reaction);
   }

   std::cerr << "Filled " << pathways_to_genes.size() << " pathways for genes.\n";
   std::cerr << "Filled " << reactions_to_genes.size() << " reactions for genes.\n";

   // Finish calculating genes to proteins
   for (const auto& entry : temp_genes_to_proteins)
      genes_to_proteins.emplace(entry.first, entry.second);

   std::unordered_set<std::string> temp_genes_mapped_to_proteins;
   std::unordered_set<std::string> temp_proteins_mapped_from_genes;
   for (const auto& entry : temp_genes_to_proteins) {
      temp_genes_mapped_to_proteins.insert(entry.first);
      temp_proteins_mapped_from_genes.insert(entry.second);
   }
   std::cerr << "Mapped " << temp_genes_mapped_to_proteins.size() << " genes to " << temp_proteins_mapped_from_genes.size() << " proteins\n";
}

void dataset::setProteinMapping(std::string_view path_file_mapping) {
   proteins = readEntities(path_file_mapping.data());
   std::cerr << "Read " << getNumProteins() << " proteins.\n";

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
   std::cerr << "Read " << getNumProteoforms() << " proteoforms.\n";

   // Calculate proteoforms to reactions and pathways
   std::ifstream map_file(path_file_mapping.data());
   std::string proteoform, reaction, pathway, leftover;

   if (!map_file.is_open())
      throw std::runtime_error("Could not open proteoform mapping file.");

   getline(map_file, leftover);  // Skip csv header line
   while (map_file.peek() != EOF) {
      proteoform = proteoform::readProteoformFromNeo4jCsv(map_file);
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

   // Calculate proteins to proteoforms
   for (const auto& proteoform : proteoforms.index_to_entities) {
      proteins_to_proteoforms.emplace(proteoform::getAccession(proteoform), proteoform);
   }
   std::unordered_set<std::string> temp_proteins_mapped_to_proteoforms;
   std::unordered_set<std::string> temp_proteoforms_mapped_from_proteins;
   for (const auto& entry : proteins_to_proteoforms) {
      temp_proteins_mapped_to_proteoforms.insert(entry.first);
      temp_proteoforms_mapped_from_proteins.insert(entry.second);
   }
   std::cerr << "Mapped " << temp_proteins_mapped_to_proteoforms.size() << " proteins to " << temp_proteoforms_mapped_from_proteins.size() << " proteoforms.\n";
}  // namespace pathway

void dataset::calculateInteractionNetworks() {
   calculateNetwork(genes, reactions_to_genes, gene_network);
   calculateNetwork(proteins, reactions_to_proteins, protein_network);
   calculateNetwork(proteoforms, reactions_to_proteoforms, proteoform_network);
}

}  // namespace pathway
