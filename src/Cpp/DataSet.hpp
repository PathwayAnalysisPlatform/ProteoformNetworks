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

#include "types.hpp"
#include "entity.hpp"
#include "proteoform.hpp"
#include "bimap.hpp"
#include "reactome.hpp"

namespace pathway {

class dataset {
  public:
   dataset() = delete;
   dataset(std::string_view path_file_gene_mapping,
           std::string_view path_file_protein_mapping,
           std::string_view path_file_proteoform_mapping);

   std::string getName() const;
   void setName(std::string_view value);

   const umss& getPathwayNames() const;

   const int getNumReactions() const;
   const int getNumPathways() const;

   const int getNumGenes() const;
   const int getNumProteins() const;
   const int getNumProteoforms() const;
   const int getNumModifiedProteins() const;
   const int getNumModifiedProteoforms() const;

   const vs& getGenes() const;
   const vs& getProteins() const;
   const vs& getProteoforms() const;
   const vs& getModifiedProteins() const;
   const vs& getModifiedProteoforms() const;

   const ummss& getGenesToReactions() const;
   const ummss& getGenesToPathways() const;
   const ummss& getProteinsToReactions() const;
   const ummss& getProteinsToPathways() const;
   const ummss& getProteoformsToReactions() const;
   const ummss& getProteoformsToPathways() const;

   const ummss& getGenesToProteins() const;
   const ummss& getProteinsToProteoforms() const;

   const ummss& getGeneNetwork() const;
   const ummss& getProteinNetwork() const;
   const ummss& getProteoformNetwork() const;

  private:
   std::string name;
   umss pathways_to_names;

   bimap phegeni_genes;
   bimap proteins;
   bimap proteoforms;
   vs modified_proteins;
   vs modified_proteoforms;

   reactome_gene_sets pathways_to_genes;
   ummss genes_to_pathways;
   um<std::string, std::bitset<REACTOME_GENES>> reactions_to_genes;
   ummss genes_to_reactions;

   um<std::string, std::bitset<REACTOME_GENES>> pathways_to_proteins;
   ummss proteins_to_pathways;
   um<std::string, std::bitset<REACTOME_GENES>> reactions_to_proteins;
   ummss proteins_to_reactions;

   um<std::string, std::bitset<REACTOME_GENES>> pathways_to_proteoforms;
   ummss proteoforms_to_pathways;
   um<std::string, std::bitset<REACTOME_GENES>> reactions_to_proteoforms;
   ummss proteoforms_to_reactions;

   ummss gene_network;
   ummss protein_network;
   ummss proteoform_network;

   ummss genes_to_proteins;
   ummss proteins_to_proteoforms;

   void setPathwayNames(std::string_view path_file_mapping);
   void setGeneMapping(std::string_view path_file_mapping);
   void setProteinMapping(std::string_view path_file_mapping);
   void setProteoformMapping(std::string_view path_file_mapping);
   void checkMappingConsistency();
   void calculateModifiedProteinsAndProteoforms();
   void calculateInteractionNetworks();

   template <size_t num_entities>
   void calculateNetwork(const bimap& entities,
                         const um<std::string, std::bitset<num_entities>>& reactions_to_entities,
                         ummss& entity_network) {
      for (const auto& reaction_entry : reactions_to_entities) {
         vs members;
         for (int I = 0; I < reaction_entry.second.size(); I++) {
            if (reaction_entry.second.test(I)) {
               members.push_back(entities.entities[I]);
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