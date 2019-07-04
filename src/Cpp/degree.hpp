#ifndef DEGREE_H_
#define DEGREE_H_

#include <fstream>
#include <string>
#include <string_view>

#include "dataset.hpp"
#include "utility.hpp"
#include "proteoform.hpp"

const std::string path_file_report_degree_analysis = "reports/degree_analysis.txt";
const std::string path_file_hits = "reports/hits.txt";
const std::string path_file_degree = "reports/degree.txt";

namespace degree {

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

void reportHits(const pathway::dataset& ds,
                std::string_view path_file_hits,
                std::string_view path_file_hits_reactions,
                std::string_view path_file_hits_pathways);

double calculateAvgDegree(entities entity_type, const vs& entities, const ummss& entity_network);

void createReportNodeDegree(const vs& entities, const ummss& entity_network, std::string_view path_file_report);

}  // namespace degree

#endif /* DEGREE_H_ */
