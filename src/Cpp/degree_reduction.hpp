#ifndef DEGREE_REDUCTION_H_
#define DEGREE_REDUCTION_H_

#include <fstream>
#include <string>
#include <string_view>

#include "dataset.hpp"

namespace degree_reduction {

struct hits_result {
   double reactions;
   double pathways;
};

void doAnalysis(const pathway::dataset& ds,
                std::string_view path_file_node_degree_genes,
                std::string_view path_file_node_degree_proteins,
                std::string_view path_file_node_degree_proteoforms);

hits_result calculateHits(pathway::entities entity_type, const pathway::dataset& ds);

double calculateAvgProteoformsPerAccession(const pathway::dataset& ds);
double calculateAvgProteoformsPerAccessionWithModification(const pathway::dataset& ds);

double calculateAvgDegree(pathway::entities entity_type, const std::vector<std::string>& entities, const std::unordered_multimap<std::string, std::string>& entity_network);

void createReportNodeDegree(const std::vector<std::string>& entities, const std::unordered_multimap<std::string, std::string>& entity_network, std::string_view path_file_report);

}  // namespace degree_reduction

#endif /* DEGREE_REDUCTION_H_ */
