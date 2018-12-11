#include "utility.hpp"

void writeFrequencies(std::string_view file_path, const std::unordered_multimap<std::string, std::string>& mapping) {
   std::cerr << "Writing " << file_path << "\n";
   std::ofstream report(file_path.data());

   if (!report.is_open()) {
      throw std::runtime_error("Problem opening frequency file.\n");
   }

   report << "OBJECT\tFREQUENCY\n";
   for (const auto& entry : mapping) {
      report << entry.first << "\t" << mapping.count(entry.first);
   }
}

void writeFrequencies(std::string_view file_path, const std::unordered_multimap<std::string, std::string>& mapping, const std::vector<std::string>& keys) {
   std::cerr << "Writing " << file_path << "\n";
   std::ofstream report(file_path.data());
   std::unordered_set<std::string> unique_keys(keys.begin(), keys.end());

   if (!report.is_open()) {
      throw std::runtime_error("Problem opening frequency file.\n");
   }

   report << "OBJECT\tFREQUENCY\n";
   for (const auto& value : keys) {
      report << value << "\t" << mapping.count(value);
   }
}
