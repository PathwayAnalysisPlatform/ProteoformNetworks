#include "degree_reduction.hpp"

namespace degree_reduction {

/* Requirements:
-- The dataset must contain the mapping for genes, proteins and proteoforms to reactions and pathways.
-- The gene mapping file should have the mapping from genes to proteins.
-- The protein mapping file should have the pathway and reaction names.*/
void doAnalysis(const pathway::dataset& ds, std::string_view report_file_path) {
   std::cout << "Degree reduction analysis...\n";

   auto avg_hits = calculateHits(pathway::entities::GENES, ds);
   std::cout << "Average reactions per gene: " << avg_hits.reactions << "\n";
   std::cout << "Average pathways per gene: " << avg_hits.pathways << "\n";

   avg_hits = calculateHits(pathway::entities::PROTEINS, ds);
   std::cout << "Average reactions per protein: " << avg_hits.reactions << "\n";
   std::cout << "Average pathways per protein: " << avg_hits.pathways << "\n";

   std::cerr << "Reactions for P31749: " << ds.getProteinsToReactions().at("P31749").size() << "\n";
   std::cerr << "Pathways for P31749: " << ds.getProteinsToPathways().at("P31749").size() << "\n";

   avg_hits = calculateHits(pathway::entities::PROTEOFORMS, ds);
   std::cout << "Average reactions per proteoform: " << avg_hits.reactions << "\n";
   std::cout << "Average pathways per proteoform: " << avg_hits.pathways << "\n";

   std::cerr << "Example accession with its proteoforms: \n";
   auto ret = ds.getProteinsToProteoforms().equal_range("P31749");
   std::cerr << "P31749 => \n";
   for (auto it = ret.first; it != ret.second; it++) {
      std::cerr << "\t" << it->second << "\n";
   }
   std::cout << "Average number of proteoforms per protein accession: " << calculateAvgProteoformsPerAccession(ds) << "\n";
   std::cout << "Average number of proteoforms per protein accession with at least a modification in any of its proteoforms: "
             << calculateAvgProteoformsPerAccessionWithModification(ds) << "\n";

}  // namespace degree_reduction

hits_result calculateHits(pathway::entities entity_type, const pathway::dataset& ds) {
   const pathway::entities_bimap* entities;
   const std::unordered_map<std::string, std::unordered_set<std::string>>* entities_to_reactions;
   const std::unordered_map<std::string, std::unordered_set<std::string>>* entities_to_pathways;

   switch (entity_type) {
      case pathway::entities::GENES:
         entities = &(ds.getGenes());
         entities_to_reactions = &(ds.getGenesToReactions());
         entities_to_pathways = &(ds.getGenesToPathways());
         break;
      case pathway::entities::PROTEINS:
         entities = &(ds.getProteins());
         entities_to_reactions = &(ds.getProteinsToReactions());
         entities_to_pathways = &(ds.getProteinsToPathways());
         break;
      case pathway::entities::PROTEOFORMS:
         entities = &(ds.getProteoforms());
         entities_to_reactions = &(ds.getProteoformsToReactions());
         entities_to_pathways = &(ds.getProteoformsToPathways());
         break;
   }

   double sum_reactions = 0.0;
   double sum_pathways = 0.0;
   for (const auto& entity : (*entities).index_to_entities) {
      sum_reactions += (*entities_to_reactions).at(entity).size();
      sum_pathways += (*entities_to_pathways).at(entity).size();
   }
   return {sum_reactions / (*entities).index_to_entities.size(), sum_pathways / (*entities).index_to_entities.size()};
}

double calculateAvgProteoformsPerAccession(const pathway::dataset& ds) {
   // Calculate average number of proteoforms for a protein
   double sum = 0.0;
   for (const auto& protein : ds.getProteins().index_to_entities) {
      sum += ds.getProteinsToProteoforms().count(protein);
   }
   return sum / ds.getNumProteins();
}

// Calculate average numer of proteoforms for proteins with at least two proteoforms with at least one modification
double calculateAvgProteoformsPerAccessionWithModification(const pathway::dataset& ds) {
   // Select accessions that have at least a proteoform with a modification
   std::unordered_set<std::string> accessions;
   for (const auto& proteoform : ds.getProteoforms().index_to_entities) {
      if (proteoform::isModified(proteoform)) {
         accessions.insert(proteoform::getAccession(proteoform));
      }
   }

   // Calculate averag proteoforms for all them
   double sum = 0.0;
   for (const auto& accession : accessions) {
      sum += ds.getProteinsToProteoforms().count(accession);
   }
   return sum / accessions.size();
}

}  // namespace degree_reduction
