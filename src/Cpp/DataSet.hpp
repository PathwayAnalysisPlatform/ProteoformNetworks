#ifndef PATHWAY_DATASET_H
#define PATHWAY_DATASET_H

#include <bitset>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "entity.hpp"
#include "proteoform.hpp"

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

   const int getNumReactions() const;
   const int getNumPathways() const;

   const int getNumGenes() const;
   const int getNumProteins() const;
   const int getNumProteoforms() const;
   const int getNumModifiedProteins() const;
   const int getNumModifiedProteoforms() const;

   const std::vector<std::string>& getGenes() const;
   const std::vector<std::string>& getProteins() const;
   const std::vector<std::string>& getProteoforms() const;
   const std::vector<std::string>& getModifiedProteins() const;
   const std::vector<std::string>& getModifiedProteoforms() const;

   const std::unordered_map<std::string, std::unordered_set<std::string>>& getGenesToReactions() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getGenesToPathways() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getProteinsToReactions() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getProteinsToPathways() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getProteoformsToReactions() const;
   const std::unordered_map<std::string, std::unordered_set<std::string>>& getProteoformsToPathways() const;

   const std::unordered_multimap<std::string, std::string>& getGenesToProteins() const;
   const std::unordered_multimap<std::string, std::string>& getProteinsToProteoforms() const;

   const std::unordered_multimap<std::string, std::string>& getGeneNetwork() const;
   const std::unordered_multimap<std::string, std::string>& getProteinNetwork() const;
   const std::unordered_multimap<std::string, std::string>& getProteoformNetwork() const;

  private:
   std::string name;
   std::unordered_map<std::string, std::string> pathways_to_names;

   entities_bimap genes;
   entities_bimap proteins;
   entities_bimap proteoforms;
   std::vector<std::string> modified_proteins;
   std::vector<std::string> modified_proteoforms;

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

   std::unordered_multimap<std::string, std::string> genes_to_proteins;
   std::unordered_multimap<std::string, std::string> proteins_to_proteoforms;

   void setPathwayNames(std::string_view path_file_mapping);
   void setGeneMapping(std::string_view path_file_mapping);
   void setProteinMapping(std::string_view path_file_mapping);
   void setProteoformMapping(std::string_view path_file_mapping);
   void checkMappingConsistency();
   void calculateModifiedProteinsAndProteoforms();
   void calculateInteractionNetworks();

   template <size_t num_entities>
   void calculateNetwork(const entities_bimap& entities,
                         const std::unordered_map<std::string, std::bitset<num_entities>>& reactions_to_entities,
                         std::unordered_multimap<std::string, std::string>& entity_network) {
      for (const auto& reaction_entry : reactions_to_entities) {
         std::vector<std::string> members;
         for (int I = 0; I < reaction_entry.second.size(); I++) {
            if (reaction_entry.second.test(I)) {
               members.push_back(entities.index_to_entities[I]);
            }
         }

         for (const auto& one_member : members) {
            for (const auto& other_member : members) {
               entity_network.emplace(one_member, other_member);
            }
         }
      }
   }
};

}  // namespace pathway

#endif /* PATHWAY_DATASET_H */