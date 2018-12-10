#include "degree_reduction.hpp"

namespace degree_reduction {

/* Requirements:
-- The dataset must contain the mapping for genes, proteins and proteoforms to reactions and pathways.
-- The gene mapping file should have the mapping from genes to proteins.
-- The protein mapping file should have the pathway and reaction names.*/
void doAnalysis(const pathway::dataset& ds,
                std::string_view path_file_node_degree_genes,
                std::string_view path_file_node_degree_proteins,
                std::string_view path_file_node_degree_proteoforms) {
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

   // Number of nodes and links in each network
   std::cout << "Gene network nodes: " << ds.getNumGenes() << " links: " << (ds.getGeneNetwork().size() / 2) << "\n";
   std::cout << "Protein network nodes: " << ds.getNumProteins() << " links: " << (ds.getProteinNetwork().size() / 2) << "\n";
   std::cout << "Proteoform network nodes: " << ds.getNumProteoforms() << " links: " << (ds.getProteoformNetwork().size() / 2) << "\n";

   std::cout << "Average degree of gene nodes: " << ds.getGeneNetwork().size() / static_cast<double>(ds.getNumGenes()) << "\n";
   std::cout << "Average degree of gene nodes: " << calculateAvgDegree(pathway::entities::GENES, ds.getGenes().index_to_entities, ds.getGeneNetwork()) << "\n";
   std::cout << "Average degree of protein nodes: " << ds.getProteinNetwork().size() / static_cast<double>(ds.getNumProteins()) << "\n";
   std::cout << "Average degree of protein nodes: " << calculateAvgDegree(pathway::entities::PROTEINS, ds.getProteins().index_to_entities, ds.getProteinNetwork()) << "\n";
   std::cout << "Average degree of proteoform nodes: " << ds.getProteoformNetwork().size() / static_cast<double>(ds.getNumProteoforms()) << "\n";
   std::cout << "Average degree of proteoform nodes: " << calculateAvgDegree(pathway::entities::PROTEOFORMS, ds.getProteoforms().index_to_entities, ds.getProteoformNetwork()) << "\n";

   // Create file with nodes and their degree. The list sorted by degree.
   createReportNodeDegree(ds.getGenes().index_to_entities, ds.getGeneNetwork(), path_file_node_degree_genes);
   createReportNodeDegree(ds.getProteins().index_to_entities, ds.getProteinNetwork(), path_file_node_degree_proteins);
   createReportNodeDegree(ds.getProteoforms().index_to_entities, ds.getProteoformNetwork(), path_file_node_degree_proteoforms);

   // Check that hub nodes reduced size.

   // Average distance L between any pair of nodes

   // Clustering coefficient

}  // namespace degree_reduction

void createReportNodeDegree(const std::vector<std::string>& entities, const std::unordered_multimap<std::string, std::string>& entity_network, std::string_view path_file_report) {
   std::ofstream file_node_degree_genes(path_file_report.data());

   if (!file_node_degree_genes.is_open()) {
      throw std::runtime_error("Could not write to node degree report file.");
   }

   std::multimap<int, std::string> freq;
   for (const auto& entity : entities) {
      freq.emplace(entity_network.count(entity), entity);
   }

   file_node_degree_genes << "ENTITY\tDEGREE\n";
   for(const auto& entry : freq){
      file_node_degree_genes << entry.second << "\t" << entry.first << "\n";
   }
}

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

double calculateAvgDegree(pathway::entities entity_type, const std::vector<std::string>& entities, const std::unordered_multimap<std::string, std::string>& entity_network) {
   double sum = 0.0;
   // For each entity of the requested type, sum its degreee
   for (const auto& entity : entities) {
      sum += entity_network.count(entity);
   }
   return sum / entities.size();
}

}  // namespace degree_reduction
