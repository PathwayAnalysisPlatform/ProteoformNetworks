#ifndef DEGREE_H_
#define DEGREE_H_

#include <fstream>
#include <string>
#include <string_view>

#include "dataset.hpp"
#include "utility.hpp"

namespace degree {

struct measures_result {
   double min;
   double max;
   double avg;
};

struct hits_result {
   double reactions;
   double pathways;
};

void doAnalysis(const pathway::dataset& ds,
                std::string_view path_file_report_degree_analysis,
                std::string_view path_file_entities,
                std::string_view path_file_proteins_per_gene,
                std::string_view path_file_proteoforms_per_protein,
                std::string_view path_file_modified_proteoforms_per_protein,
                std::string_view path_file_hits,
                std::string_view path_file_degree);

void reportEntities(const pathway::dataset& ds,
                    std::string_view path_file_entities,
                    std::string_view path_file_proteins_per_gene,
                    std::string_view path_file_proteoforms_per_protein,
                    std::string_view path_file_modified_proteoforms_per_protein);

void reportHits(std::string_view path_file_report_degree_analysis, std::string_view path_file_hits);
const measures_result calculateMeasures(const std::unordered_multimap<std::string, std::string>& mapping);
const measures_result calculateMeasuresWithSelectedKeys(const std::unordered_multimap<std::string, std::string>& mapping,
                                                        const std::vector<std::string>& keys);

const hits_result calculateHits(pathway::entities entity_type, const pathway::dataset& ds);

double calculateAvgDegree(pathway::entities entity_type, const std::vector<std::string>& entities, const std::unordered_multimap<std::string, std::string>& entity_network);

void createReportNodeDegree(const std::vector<std::string>& entities, const std::unordered_multimap<std::string, std::string>& entity_network, std::string_view path_file_report);

}  // namespace degree

#endif /* DEGREE_H_ */
