#ifndef PATHWAY_DATASET_H
#define PATHWAY_DATASET_H

#include <bitset>
#include <fstream>
#include <string>
#include <string_view>
#include <unordered_map>

#include "entity.h"

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

   void setGeneSets(std::string_view path_file_mapping);
   void setProteinSets(std::string_view path_file_mapping);
   void setProteoformSets(std::string_view path_file_mapping);
   void setPathwayNames(std::string_view path_file_mapping);
   void calculateInteractionNetworks();

  private:
   std::string name;
   std::unordered_map<std::string, std::string> pathways_to_names;

   entities_bimap genes;
   entities_bimap proteins;
   entities_bimap proteoforms;

   std::unordered_map<std::string, std::bitset<NUM_GENES>> gene_mapping;
   std::unordered_map<std::string, std::bitset<NUM_PROTEINS>> protein_mapping;
   std::unordered_map<std::string, std::bitset<NUM_PROTEOFORMS>> proteoform_mapping;

   std::unordered_multimap<std::string, std::string> gene_network;
   std::unordered_multimap<std::string, std::string> protein_network;
   std::unordered_multimap<std::string, std::string> proteoform_network;
};

}  // namespace pathway

#endif /* PATHWAY_DATASET_H */