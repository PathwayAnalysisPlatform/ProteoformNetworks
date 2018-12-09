#ifndef DEGREE_REDUCTION_H_
#define DEGREE_REDUCTION_H_

#include <string>
#include <string_view>
#include "dataset.hpp"

namespace degree_reduction {

struct hits_result {
    double reactions;
    double pathways;
};

void doAnalysis(const pathway::dataset& ds, std::string_view report_file_path);

hits_result calculateHits(pathway::entities entity_type, const pathway::dataset& ds);

double calculateAvgProteoformsPerAccession(const pathway::dataset& ds);
double calculateAvgProteoformsPerAccessionWithModification(const pathway::dataset& ds);

}

#endif /* DEGREE_REDUCTION_H_ */
