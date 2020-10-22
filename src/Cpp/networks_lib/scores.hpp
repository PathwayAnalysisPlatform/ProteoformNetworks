#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#define PHEGENI_HPP

#include <assert.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <string_view>
#include <bitset.h>
#include <bitset>
#include <unordered_map>
#include <utility>

#include "types.hpp"
#include "bimap_str_int.hpp"
#include "overlap_types.hpp"

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

// Calculate score between al pairs of bitsets
// The sets are the second value of each entry in the sets parameter.
// The score is a function capable of calculating the overlap with bitsets.
// Returns only the sets within the module sizes and with a score greater than 0.
pair_map<double>
getScores(const vb &vertex_sets, std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>)> score_function,
          const int min_module_size, const int max_module_size);

// Calculate score between the selected pairs.
// The sets are the second value of each entry in the sets parameter.
// The score is a function capable of calculating the overlap with bitsets.
pair_map<double>
getScores(const vb &vertex_sets, std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>)> score_function,
          const pair_map<double> &prev_score);

pair_map<double>
getScores(const vb &vertex_sets,
          const vusi &edges,
          std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>, vusi)> score_function,
          const pair_map<double> &prev_score);

double calculate_interface_size_nodes(const base::dynamic_bitset<> &V1,
                                      const base::dynamic_bitset<> &V2,
                                      const vusi &E);

double calculate_interface_size_edges(const base::dynamic_bitset<> &V1,
                                      const base::dynamic_bitset<> &V2,
                                      const vusi &E);

#endif /* UTILITY_HPP_ */
