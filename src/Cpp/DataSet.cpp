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
   calculateModifiedProteinsAndProteoforms();

   // Check consistency of reactions and pathways for the three types of entities
   checkMappingConsistency();

   // Create interaction network
   calculateInteractionNetworks();
}

std::string dataset::getName() const {
   return name;
}

void dataset::setName(std::string_view value) {
   name = value;
}

const umss& dataset::getPathwayNames() const {
   return pathways_to_names;
}

const int dataset::getNumReactions() const {
   return reactions_to_genes.size();
}

const int dataset::getNumPathways() const {
   return pathways_to_genes.size();
}

const int dataset::getNumGenes() const {
   return phegeni_genes.entities.size();
}
const int dataset::getNumProteins() const {
   return proteins.entities.size();
}
const int dataset::getNumProteoforms() const {
   return proteoforms.entities.size();
}

const int dataset::getNumModifiedProteins() const {
   return modified_proteins.size();
}
const int dataset::getNumModifiedProteoforms() const {
   return modified_proteoforms.size();
}

const vs& dataset::getGenes() const {
   return phegeni_genes.entities;
}
const vs& dataset::getProteins() const {
   return proteins.entities;
}
const vs& dataset::getProteoforms() const {
   return proteoforms.entities;
}
const vs& dataset::getModifiedProteins() const {
   return modified_proteins;
}
const vs& dataset::getModifiedProteoforms() const {
   return modified_proteoforms;
}

const ummss& dataset::getGenesToReactions() const {
   return genes_to_reactions;
}
const ummss& dataset::getGenesToPathways() const {
   return genes_to_pathways;
}
const ummss& dataset::getProteinsToReactions() const {
   return proteins_to_reactions;
}
const ummss& dataset::getProteinsToPathways() const {
   return proteins_to_pathways;
}
const ummss& dataset::getProteoformsToReactions() const {
   return proteoforms_to_reactions;
}
const ummss& dataset::getProteoformsToPathways() const {
   return proteoforms_to_pathways;
}
const ummss& dataset::getGenesToProteins() const {
   return genes_to_proteins;
}
const ummss& dataset::getProteinsToProteoforms() const {
   return proteins_to_proteoforms;
}
const ummss& dataset::getGeneNetwork() const {
   return gene_network;
}
const ummss& dataset::getProteinNetwork() const {
   return protein_network;
}
const ummss& dataset::getProteoformNetwork() const {
   return proteoform_network;
}

void dataset::setPathwayNames(std::string_view path_file_mapping) {
   std::cerr << "Loading pathways\n";
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
}

void dataset::setGeneMapping(std::string_view path_file_mapping) {
   std::cerr << "Loading gene mapping\n";
   phegeni_genes = createBimap(path_file_mapping);
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

      pathways_to_genes[pathway].set(phegeni_genes.indexes.at(gene));
      genes_to_pathways.emplace(gene, pathway);

      reactions_to_genes[reaction].set(phegeni_genes.indexes.at(gene));
      genes_to_reactions.emplace(gene, reaction);
   }

   // Finish calculating genes to proteins
   for (const auto& entry : temp_genes_to_proteins)
      genes_to_proteins.emplace(entry.first, entry.second);
}

void dataset::setProteinMapping(std::string_view path_file_mapping) {
   std::cerr << "Loading protein mapping\n";
   proteins = createBimap(path_file_mapping.data());
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

      pathways_to_proteins[pathway].set(proteins.indexes.at(protein));
      proteins_to_pathways.emplace(protein, pathway);

      reactions_to_proteins[reaction].set(proteins.indexes.at(protein));
      proteins_to_reactions.emplace(protein, reaction);
   }
}

void dataset::setProteoformMapping(std::string_view path_file_mapping) {
   std::cerr << "Loading proteoform mapping\n";
   proteoforms = createBimap(path_file_mapping);
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

      pathways_to_proteoforms[pathway].set(proteoforms.indexes.at(proteoform));
      proteoforms_to_pathways.emplace(proteoform, pathway);

      reactions_to_proteoforms[reaction].set(proteoforms.indexes.at(proteoform));
      proteoforms_to_reactions.emplace(proteoform, reaction);
   }

   // Calculate proteins to proteoforms
   for (const auto& proteoform : proteoforms.entities) {
      proteins_to_proteoforms.emplace(proteoform::getAccession(proteoform), proteoform);
   }
}  // namespace pathway

void dataset::calculateModifiedProteinsAndProteoforms() {
   uss temp_modified_proteins;
   for (const auto& proteoform : getProteoforms()) {
      if (proteoform::isModified(proteoform)) {
         temp_modified_proteins.insert(proteoform::getAccession(proteoform));
         modified_proteoforms.push_back(proteoform);
      }
   }
   for (const auto& accession : temp_modified_proteins) {
      modified_proteins.push_back(accession);
   }
}

void dataset::calculateInteractionNetworks() {
   std::cerr << "Calculating networks...\n";
   calculateNetwork(phegeni_genes, reactions_to_genes, gene_network);
   calculateNetwork(proteins, reactions_to_proteins, protein_network);
   calculateNetwork(proteoforms, reactions_to_proteoforms, proteoform_network);
}

void dataset::checkMappingConsistency() {
   std::cerr << "Filled " << pathways_to_genes.size() << " pathways for genes.\n";
   std::cerr << "Filled " << reactions_to_genes.size() << " reactions for genes.\n";

   std::cerr << "Filled " << pathways_to_proteins.size() << " pathways for proteins.\n";
   std::cerr << "Filled " << reactions_to_proteins.size() << " reactions for proteins.\n";

   std::cerr << "Filled " << pathways_to_proteoforms.size() << " pathways for proteoforms.\n";
   std::cerr << "Filled " << reactions_to_proteoforms.size() << " reactions for proteoforms.\n";

   if (reactions_to_genes.size() != reactions_to_proteins.size() != reactions_to_proteoforms.size()) {
      throw std::runtime_error("The number of reactions mapped for each entity type differs.\n");
   }

   if (pathways_to_genes.size() != pathways_to_proteins.size() != pathways_to_proteoforms.size()) {
      throw std::runtime_error("The number of pathways mapped for each entity type differs.\n");
   }

   uss temp_genes_mapped_to_proteins;
   uss temp_proteins_mapped_from_genes;
   for (const auto& entry : genes_to_proteins) {
      temp_genes_mapped_to_proteins.insert(entry.first);
      temp_proteins_mapped_from_genes.insert(entry.second);
   }
   std::cerr << "Mapped " << temp_genes_mapped_to_proteins.size() << " genes to " << temp_proteins_mapped_from_genes.size() << " proteins\n";

   if (getNumGenes() != temp_genes_mapped_to_proteins.size()) {
      throw std::runtime_error("The number of genes mapped to protein differs with the total number of genes.\n");
   }

   if (getNumProteins() != temp_proteins_mapped_from_genes.size()) {
      throw std::runtime_error("The number of proteins mapped from genes differs with the total number of proteins.");
   }

   uss temp_proteins_mapped_to_proteoforms;
   uss temp_proteoforms_mapped_from_proteins;
   for (const auto& entry : proteins_to_proteoforms) {
      temp_proteins_mapped_to_proteoforms.insert(entry.first);
      temp_proteoforms_mapped_from_proteins.insert(entry.second);
   }
   std::cerr << "Mapped " << temp_proteins_mapped_to_proteoforms.size() << " proteins to " << temp_proteoforms_mapped_from_proteins.size() << " proteoforms.\n";

   if (getNumProteins() != temp_proteins_mapped_to_proteoforms.size()) {
      throw std::runtime_error("The number of proteins mapped to proteoforms differs from the total number of proteins.\n");
   }

   if (getNumProteoforms() != temp_proteoforms_mapped_from_proteins.size()) {
      throw std::runtime_error("The number of proteoforms mapped from proteins differs from the total number of proteoforms.\n");
   }
}

}  // namespace pathway
