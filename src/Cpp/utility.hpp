#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

void writeFrequencies(std::string_view file_path, const std::unordered_multimap<std::string, std::string>& mapping);
void writeFrequencies(std::string_view file_path, const std::unordered_multimap<std::string, std::string>& mapping, const std::vector<std::string>& keys);

#endif /* UTILITY_HPP */