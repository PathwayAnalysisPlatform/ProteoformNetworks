#ifndef MODIFIED_OVERLAP_H_
#define MODIFIED_OVERLAP_H_

#include "../overlap.hpp"

#include <regex>

const float MIN_MODIFIED_ALL_MEMBERS_RATIO = 0.1;
const float MIN_MODIFIED_OVERLAP_MEMBERS_RATIO = 0.9;

const std::regex RGX_MODIFICATION{"[;,]\\d{5}"};

void doModifiedOverlapAnalysis(const std::string& path_file_proteoform_search,
                               const std::string& path_file_report,
                               const std::string& path_file_proteins,
                               const std::string& path_file_proteoforms,
                               const std::string& path_file_modifications);

#endif /* MODIFIED_OVERLAP_H_ */