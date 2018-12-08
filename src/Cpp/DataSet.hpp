#ifndef PATHWAY_DATASET_H
#define PATHWAY_DATASET_H

#include <bitset>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include "entity.hpp"

namespace pathway {

const size_t NUM_GENES = 23970;
const size_t NUM_PROTEINS = 10778;
const size_t NUM_PROTEOFORMS = 13911;

class dataset {
  public:
   dataset() = delete;
   dataset(std::string_view path_file_gene_mapping,
           std::string_view path_file_protein_mapping,
           std::string_view path_file_proteoform_mapping);

   std::string getName() const;
   void setName(std::string_view value);

   const std::unordered_map<std::string, std::string>& getPathwayNames() const;

   const entities_bimap& getGenes() const;
   const entities_bimap& getProteins() const;
   const entities_bimap& getProteoforms() const;

   const std::unordered_map<std::string, std::unordered_set<std::string>>& getGenesToReactions() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getGenesToPathways() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getProteinsToReactions() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getProteinsToPathways() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getProteoformsToReactions() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getProteoformsToPathways() const;

  private:
   std::string name;
   std::unordered_map<std::string, std::string> pathways_to_names;

   entities_bimap genes;
   entities_bimap proteins;
   entities_bimap proteoforms;

   std::unordered_map<std::string, std::bitset<NUM_GENES>> pathways_to_genes;
   std::unordered_map<std::string, std::unordered_set<std::string>> genes_to_pathways;
   std::unordered_map<std::string, std::bitset<NUM_GENES>> reactions_to_genes;
   std::unordered_map<std::string, std::unordered_set<std::string>> genes_to_reactions;

   std::unordered_map<std::string, std::bitset<NUM_GENES>> pathways_to_proteins;
   std::unordered_map<std::string, std::unordered_set<std::string>> proteins_to_pathways;
   std::unordered_map<std::string, std::bitset<NUM_GENES>> reactions_to_proteins;
   std::unordered_map<std::string, std::unordered_set<std::string>> proteins_to_reactions;

   std::unordered_map<std::string, std::bitset<NUM_GENES>> pathways_to_proteoforms;
   std::unordered_map<std::string, std::unordered_set<std::string>> proteoforms_to_pathways;
   std::unordered_map<std::string, std::bitset<NUM_GENES>> reactions_to_proteoforms;
   std::unordered_map<std::string, std::unordered_set<std::string>> proteoforms_to_reactions;

   std::unordered_multimap<std::string, std::string> gene_network;
   std::unordered_multimap<std::string, std::string> protein_network;
   std::unordered_multimap<std::string, std::string> proteoform_network;

   void setPathwayNames(std::string_view path_file_mapping);
   void setGeneMapping(std::string_view path_file_mapping);
   void setProteinMapping(std::string_view path_file_mapping);
   void setProteoformMapping(std::string_view path_file_mapping);
   void calculateInteractionNetworks();
};

}  // namespace pathway

#endif /* PATHWAY_DATASET_H */