#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#define PHEGENI_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <string_view>
#include <bitset.h>
#include <bitset>
#include <unordered_map>
#include <utility>

#include "types.hpp"
#include "phegeni.hpp"
#include "bimap_str_int.hpp"

struct measures_result {
    double min;
    double max;
    double avg;
};

const measures_result calculateMeasures(const ummss &mapping);

const measures_result calculateMeasuresWithSelectedKeys(const ummss &mapping, const vs &keys);

void writeFrequencies(std::string_view file_path, const ummss &mapping);

void writeFrequencies(std::string_view file_path, const ummss &mapping, const vs &keys);

void
writeMeasures(std::ofstream &report, const measures_result &measures, std::string_view label1, std::string_view label2);

void writeMeasures(std::ofstream &report, const ummss &mapping, std::string_view label1, std::string_view label2);

// Calculate Jaccard index, which is intersection over union
double getJaccardSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

double getOverlapSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

double getOverlapSize(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

// Calculate similarity score between al pairs of bitsets
// The sets are the second value of each entry in the sets parameter.
// The score is a function capable of calculating the overlap with bitsets.
pair_map<double>
getScores(const vb &sets, std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>)> score_function);

void writeScores(const bimap_str_int &groups, const modules &entity_modules,
                 const pair_map<double> &scores, const pair_map<double> &overlap_sizes, std::string_view path_output);

#endif /* UTILITY_HPP_ */
