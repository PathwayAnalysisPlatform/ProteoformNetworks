#include "utility.hpp"

vs convert_uss_to_vs(const uss& a_set) {
	vs result;
	result.assign(a_set.begin(), a_set.end());
	sort(result.begin(), result.end());
	return result;
}

const measures_result calculateMeasures(const ummss& mapping) {
   double min = static_cast<double>(mapping.count(mapping.begin()->first));
   double max = static_cast<double>(mapping.count(mapping.begin()->first));
   double avg = 0.0;
   double sum = 0.0;
   for (const auto& entry : mapping) {
      double freq = static_cast<double>(mapping.count(entry.first));
      if (freq < min)
         min = freq;
      if (freq > max)
         max = freq;
      sum += freq;
   }
   avg = static_cast<double>(sum) / mapping.size();
   return {min, max, avg};
}

// Calculate average numer of proteoforms for proteins with at least two proteoforms with at least one modification
const measures_result calculateMeasuresWithSelectedKeys(const ummss& mapping,
                                                        const vs& keys) {
   // Calculate averag proteoforms for all them
   double min = static_cast<double>(mapping.count(*keys.begin()));
   double max = static_cast<double>(mapping.count(*keys.begin()));
   double avg = 0.0;
   double sum = 0.0;
   for (const auto& value : keys) {
      double freq = static_cast<double>(mapping.count(value));
      if (freq < min)
         min = freq;
      if (freq > max)
         max = freq;
      sum += mapping.count(value);
   }
   avg = sum / keys.size();
   return {min, max, avg};
}

void writeFrequencies(std::ofstream& report, const ummss& mapping, std::string_view label) {
   if (!report.is_open()) {
      throw std::runtime_error("Problem opening frequency file.\n");
   }

   for (const auto& entry : mapping) {
      report << entry.first << "\t" << mapping.count(entry.first);
   }
}

void writeFrequencies(std::string_view file_path, const ummss& mapping) {
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

void writeFrequencies(std::string_view file_path, const ummss& mapping, const vs& keys) {
   std::cerr << "Writing " << file_path << "\n";
   std::ofstream report(file_path.data());
   uss unique_keys(keys.begin(), keys.end());

   if (!report.is_open()) {
      throw std::runtime_error("Problem opening frequency file.\n");
   }

   report << "OBJECT\tFREQUENCY\n";
   for (const auto& value : keys) {
      report << value << "\t" << mapping.count(value);
   }
}


void writeMeasures(std::ofstream& report, const measures_result& measures, std::string_view label1, std::string_view label2) {
   if (!report.is_open()) {
      throw std::runtime_error("Could not write measures to report.");
   }
   report << "Average " << label1 << " per " << label2 << ": " << measures.avg << "\n";
   report << "Min: " << measures.min << "\n";
   report << "Max: " << measures.max << "\n";
}

void writeMeasures(std::ofstream& report, const ummss& mapping, std::string_view label1, std::string_view label2) {
   if (!report.is_open()) {
      throw std::runtime_error("Could not write measures to report.");
   }
   auto measures = calculateMeasures(mapping);
   writeMeasures(report, measures, label1, label2);
}