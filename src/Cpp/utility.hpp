#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>

#include "types.hpp"

vs convert_uss_to_vs(const uss& a_set);

struct measures_result {
	double min;
	double max;
	double avg;
};

const measures_result calculateMeasures(const ummss& mapping);
const measures_result calculateMeasuresWithSelectedKeys(const ummss& mapping,
	const vs& keys);

void writeFrequencies(std::string_view file_path, const ummss& mapping);
void writeFrequencies(std::string_view file_path, const ummss& mapping, const vs& keys);

void writeMeasures(std::ofstream& report, const measures_result& measures, std::string_view label1, std::string_view label2);
void writeMeasures(std::ofstream& report, const ummss& mapping, std::string_view label1, std::string_view label2);

#endif /* UTILITY_HPP_ */