#ifndef DEGREE_REDUCTION_H_
#define DEGREE_REDUCTION_H_

#include <string>
#include <string_view>
#include "dataset.hpp"

namespace degree_reduction {

void doAnalysis(pathway::dataset pathwayDataSet, std::string_view report_file_path);

}

#endif /* DEGREE_REDUCTION_H_ */
