#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <string_view>
#include <bitset.h>

#include "../others/types.hpp"
#include "phegeni.hpp"

struct measures_result {
    double min;
    double max;
    double avg;
};

const measures_result calculateMeasures(const ummss& mapping);
const measures_result calculateMeasuresWithSelectedKeys(const ummss& mapping, const vs& keys);

void writeFrequencies(std::string_view file_path, const ummss& mapping);
void writeFrequencies(std::string_view file_path, const ummss& mapping, const vs& keys);

void writeMeasures(std::ofstream& report, const measures_result& measures, std::string_view label1, std::string_view label2);
void writeMeasures(std::ofstream& report, const ummss& mapping, std::string_view label1, std::string_view label2);

// Calculate Jaccard index, which is intersection over union
double getJaccardSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

double getOverlapSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

#endif /* UTILITY_HPP_ */
