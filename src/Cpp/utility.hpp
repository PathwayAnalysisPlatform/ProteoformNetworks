#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template<typename K, typename V = K>
using mapping = std::multimap<K, V>;
using string_mapping = mapping<std::string>;

struct measures_result {
   double min;
   double max;
   double avg;
};

const measures_result calculateMeasures(const std::unordered_multimap<std::string, std::string>& mapping);
const measures_result calculateMeasuresWithSelectedKeys(const std::unordered_multimap<std::string, std::string>& mapping,
                                                        const std::vector<std::string>& keys);

void writeFrequencies(std::string_view file_path, const std::unordered_multimap<std::string, std::string>& mapping);
void writeFrequencies(std::string_view file_path, const std::unordered_multimap<std::string, std::string>& mapping, const std::vector<std::string>& keys);

void writeMeasures(std::ofstream& report, const measures_result& measures, std::string_view label1, std::string_view label2);
void writeMeasures(std::ofstream& report, const std::unordered_multimap<std::string, std::string>& mapping, std::string_view label1, std::string_view label2);

#endif /* UTILITY_HPP */