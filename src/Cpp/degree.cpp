#include "degree.hpp"

namespace degree {

/* Requirements:
-- The dataset must contain the mapping for genes, proteins and proteoforms to reactions and pathways.
-- The gene mapping file should have the mapping from genes to proteins.
-- The protein mapping file should have the pathway and reaction names.*/
void doAnalysis(const pathway::dataset& ds,
                std::string_view path_file_report_degree_analysis,
                std::string_view path_file_entities,
                std::string_view path_file_proteins_per_gene,
                std::string_view path_file_proteoforms_per_protein,
                std::string_view path_file_modified_proteoforms_per_protein,
                std::string_view path_file_node_degree_genes,
                std::string_view path_file_node_degree_proteins,
                std::string_view path_file_node_degree_proteoforms,
                std::string_view path_file_hits,
                std::string_view path_file_hits_reactions,
                std::string_view path_file_hits_pathways,
                std::string_view path_file_degree) {
   std::cout << "Degree analysis...\n";
   reportEntities(ds, path_file_entities,
                  path_file_proteins_per_gene,
                  path_file_proteoforms_per_protein,
                  path_file_modified_proteoforms_per_protein);

   reportHits(ds, path_file_report_degree_analysis, path_file_hits, path_file_hits_reactions);

   // Number of nodes and links in each network
   std::cout << "Gene network nodes: " << ds.getNumGenes() << " links: " << (ds.getGeneNetwork().size() / 2) << "\n";
   std::cout << "Protein network nodes: " << ds.getNumProteins() << " links: " << (ds.getProteinNetwork().size() / 2) << "\n";
   std::cout << "Proteoform network nodes: " << ds.getNumProteoforms() << " links: " << (ds.getProteoformNetwork().size() / 2) << "\n";

   std::cout << "Average degree of gene nodes: " << ds.getGeneNetwork().size() / static_cast<double>(ds.getNumGenes()) << "\n";
   std::cout << "Average degree of protein nodes: " << ds.getProteinNetwork().size() / static_cast<double>(ds.getNumProteins()) << "\n";
   std::cout << "Average degree of proteoform nodes: " << ds.getProteoformNetwork().size() / static_cast<double>(ds.getNumProteoforms()) << "\n";

   // Create file with nodes and their degree. The list sorted by degree.
   writeFrequencies(path_file_node_degree_genes, ds.getGeneNetwork());
   writeFrequencies(path_file_node_degree_proteins, ds.getProteinNetwork());
   writeFrequencies(path_file_node_degree_proteoforms, ds.getProteoformNetwork());

   // Check which hub nodes reduced size
   // Report node degree
   std::ofstream file_report_degree("reports/degree.txt");

   if (!file_report_degree.is_open()) {
      throw std::runtime_error("Could not open degree report file.\n");
   }

   file_report_degree << "GENE\tGENE_DEGREE\tPROTEIN\tPROTEIN_DEGREE\tPROTEOFORM\tPROTEOFORM_DEGREE\tVARIATION_GENE_TO_PROTEIN\tVARIATION_PROTEIN_TO_PROTEOFORM\n";
   for (const auto& gene_entry : ds.getGenesToProteins()) {  // For each gene
      auto ret = ds.getProteinsToProteoforms().equal_range(gene_entry.second);
      for (auto it = ret.first; it != ret.second; it++) {  // For each protein product of that gene
         std::string_view gene = gene_entry.first;
         std::string_view protein = gene_entry.second;
         std::string_view proteoform = it->second;
         size_t gene_degree = ds.getGeneNetwork().count(gene.data());
         size_t protein_degree = ds.getProteinNetwork().count(protein.data());
         size_t proteoform_degree = ds.getProteoformNetwork().count(proteoform.data());
         double variation_gene_to_protein = (protein_degree / gene_degree) - 1.0;
         double variation_protein_to_proteoform = (proteoform_degree / protein_degree) - 1.0;

         file_report_degree << gene << "\t" << gene_degree << "\t";
         file_report_degree << protein << "\t" << protein_degree << "\t";
         file_report_degree << proteoform << "\t" << proteoform_degree << "\t";
         file_report_degree << variation_gene_to_protein << "\t" << variation_protein_to_proteoform << "\n";

         if (variation_gene_to_protein > 0) {
            throw std::runtime_error("Protein node with higher degree as the respective gene node degree.");
         }
      }
   }

   std::cerr << "Example accession with its proteoforms: \n";
   auto ret = ds.getProteinsToProteoforms().equal_range("P31749");
   std::cerr << "P31749 => " << ds.getProteinNetwork().count("P31749") << "\n";
   for (auto it = ret.first; it != ret.second; it++) {
      std::cerr << "\t" << it->second << " => " << ds.getProteoformNetwork().count(it->second) << "\n";
   }

   // Clustering coefficient
   // Calculate clustering for each pathway and check which varied the most

   // Average distance L between any pair of nodes
   // Calculate the average distance between nodes in each pathway. Check for pathways with more variation

}  // namespace degree

void reportEntities(const pathway::dataset& ds,
                    std::string_view path_file_entities,
                    std::string_view path_file_proteins_per_gene,
                    std::string_view path_file_proteoforms_per_protein,
                    std::string_view path_file_modified_proteoforms_per_protein) {
   std::ofstream report_entities(path_file_entities.data());
   std::ofstream report_proteins_per_gene(path_file_proteins_per_gene.data());
   std::ofstream report_proteoforms_per_protein(path_file_proteoforms_per_protein.data());

   if (!report_entities.is_open()) {
      throw std::runtime_error("Could not open report for entities file.\n");
   }
   if (!report_proteins_per_gene.is_open()) {
      throw std::runtime_error("Could not open report for proteins per gene.\n");
   }
   if (!report_proteoforms_per_protein.is_open()) {
      throw std::runtime_error("Could not open report for proteoforms per protein.\n");
   }

   // Number of reactions and pathways
   report_entities << "Number of reactions: " << ds.getNumReactions() << "\n";
   report_entities << "Number of pathways: " << ds.getNumPathways() << "\n";

   // Number of genes, proteins and proteoforms
   report_entities << "Number of genes: " << ds.getNumGenes() << "\n";
   report_entities << "Number of proteins: " << ds.getNumProteins() << "\n";
   report_entities << "Number of proteoforms: " << ds.getNumProteoforms() << "\n";

   // Number of proteins per gene: average, max, min
   auto measures = calculateMeasures(ds.getGenesToProteins());
   report_entities << "\n";
   if (measures.avg != static_cast<double>(ds.getNumProteins()) / ds.getNumGenes())
      throw std::runtime_error("The average number of accessions per gene has a problem with the calculation.\n");
   writeMeasures(report_entities, measures, "proteins", "gene");
   writeFrequencies(path_file_proteins_per_gene, ds.getGenesToProteins());

   // Number of proteoforms per protein: average, max, min
   measures = calculateMeasures(ds.getProteinsToProteoforms());
   report_entities << "\n";
   if (measures.avg != static_cast<double>(ds.getNumProteoforms()) / ds.getNumProteins())
      throw std::runtime_error("There is a problem in the calculation of proteoforms per accession.\n");
   writeMeasures(report_entities, measures, "proteoforms", "protein");
   writeFrequencies(path_file_proteoforms_per_protein, ds.getProteinsToProteoforms());

   measures = calculateMeasuresWithSelectedKeys(ds.getProteinsToProteoforms(), ds.getModifiedProteins());
   writeMeasures(report_entities, measures, "proteoforms", "modified protein");
   writeFrequencies(path_file_modified_proteoforms_per_protein, ds.getProteinsToProteoforms(), ds.getModifiedProteins());

   // Number of modifications per proteoform: average, max, min
   measures = proteoform::calculateModificationsPerProteoform(ds.getProteoforms());
   writeMeasures(report_entities, measures, "modifications", "proteoform");

   // Frequency of modifications
   //TODO

   // Frequency of modification types
   //TODO
}

void reportHits(const pathway::dataset& ds,
                std::string_view path_file_hits,
                std::string_view path_file_hits_reactions,
                std::string_view path_file_hits_pathways) {
   std::ofstream report_hits(path_file_hits.data());
   std::ofstream report_hits_reactions(path_file_hits_reactions.data());
   std::ofstream report_hits_pathways(path_file_hits_pathways.data());

   writeMeasures(report_hits, ds.getGenesToReactions(), "reactions", "gene");
   //writeFrequencies(path_file);
   writeMeasures(report_hits, ds.getGenesToPathways(), "pathways", "gene");
   writeMeasures(report_hits, ds.getProteinsToReactions(), "reactions", "protein");
   writeMeasures(report_hits, ds.getProteoformsToReactions(), "reactions", "proteoform");
   writeMeasures(report_hits, ds.getProteoformsToPathways(), "pathways", "proteoform");

   std::cerr << "Reactions for P31749: " << ds.getProteinsToReactions().count("P31749") << "\n";
   std::cerr << "Pathways for P31749: " << ds.getProteinsToPathways().count("P31749") << "\n";
}

double calculateAvgDegree(entities entity_type, const vs& entities, const ummss& entity_network) {
   double sum = 0.0;
   // For each entity of the requested type, sum its degreee
   for (const auto& entity : entities) {
      sum += entity_network.count(entity);
   }
   return sum / entities.size();
}

}  // namespace degree
