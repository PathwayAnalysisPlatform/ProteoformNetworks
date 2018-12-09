#include "degree_reduction.hpp"

using namespace std;

namespace degree_reduction {

/* Requirements:
-- The dataset must contain the mapping for genes, proteins and proteoforms to reactions and pathways.
-- The gene mapping file should have the mapping from genes to proteins.
-- The protein mapping file should have the pathway and reaction names.*/
void doAnalysis(pathway::dataset dataset_reactome, string_view report_file_path) {
   std::cout << "Degree reduction analysis...\n";

   // Calculate average number of reactions and pathways where a gene, protein and proteoform participate
   double avg_reactions_genes;
   double avg_pathways_genes;
   double avg_reactions_proteins;
   double avg_pathways_proteins;
   double avg_reactions_proteoforms;
   double avg_pathways_proteoforms;

   double sum_reactions = 0.0;
   double sum_pathways = 0.0;
   for (const auto& gene : dataset_reactome.getGenes().index_to_entities) {
      sum_reactions += dataset_reactome.getGenesToReactions().at(gene).size();
      sum_pathways += dataset_reactome.getGenesToPathways().at(gene).size();
   }
   avg_reactions_genes = sum_reactions / dataset_reactome.getNumGenes();
   avg_pathways_genes = sum_pathways / dataset_reactome.getNumGenes();
   std::cerr << "Average reactions per gene: " << avg_reactions_genes << "\n";
   std::cerr << "Average pathways per gene: " << avg_pathways_genes << "\n";

   sum_reactions = 0.0;
   sum_pathways = 0.0;
   for (const auto& protein : dataset_reactome.getProteins().index_to_entities) {
      sum_reactions += dataset_reactome.getProteinsToReactions().at(protein).size();
      sum_pathways += dataset_reactome.getProteinsToPathways().at(protein).size();
   }
   avg_reactions_proteins = sum_reactions / dataset_reactome.getNumProteins();
   avg_pathways_proteins = sum_pathways / dataset_reactome.getNumProteins();
   std::cerr << "Average reactions per protein: " << avg_reactions_proteins << "\n";
   std::cerr << "Average pathways per protein: " << avg_pathways_proteins << "\n";

   std::cerr << "Reactions for P31749: " << dataset_reactome.getProteinsToReactions().at("P31749").size() << "\n";
   std::cerr << "Pathways for P31749: " << dataset_reactome.getProteinsToPathways().at("P31749").size() << "\n";

   sum_reactions = 0.0;
   sum_pathways = 0.0;
   for (const auto& proteoform : dataset_reactome.getProteoforms().index_to_entities) {
      sum_reactions += dataset_reactome.getProteoformsToReactions().at(proteoform).size();
      sum_pathways += dataset_reactome.getProteoformsToPathways().at(proteoform).size();
   }
   avg_reactions_proteoforms = sum_reactions / dataset_reactome.getNumProteoforms();
   avg_pathways_proteoforms = sum_pathways / dataset_reactome.getNumProteoforms();
   std::cerr << "Average reactions per proteoform: " << avg_reactions_proteoforms << "\n";
   std::cerr << "Average pathways per proteoform: " << avg_pathways_proteoforms << "\n";

   // Calculate average number of proteoforms for a protein
   double avg_proteoforms_per_protein;
   double sum_proteoforms = 0.0;
   for (const auto& protein : dataset_reactome.getProteins().index_to_entities) {
      sum_proteoforms += dataset_reactome.getProteinsToProteoforms().count(protein);
   }
   avg_proteoforms_per_protein = sum_proteoforms / dataset_reactome.getNumProteins();
   std::cerr << "Average number of proteoforms per protein accession: " << avg_proteoforms_per_protein << "\n";

   std::cerr << "Example accession with its proteoforms: \n";
   auto ret = dataset_reactome.getProteinsToProteoforms().equal_range("P31749");
   std::cerr << "P31749 => \n";
   for (auto it = ret.first; it != ret.second; it++) {
      cout << "\t" << it->second << "\n";
   }

}  // namespace degree_reduction

}  // namespace degree_reduction
