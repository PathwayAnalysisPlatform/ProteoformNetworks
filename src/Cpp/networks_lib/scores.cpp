#include "scores.hpp"

using namespace std;

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

void writeFrequencies(ofstream& report, const ummss& mapping, string_view label) {
    if (!report.is_open()) {
        throw runtime_error("Problem opening frequency file.\n");
    }

    for (const auto& entry : mapping) {
        report << entry.first << "\t" << mapping.count(entry.first);
    }
}

void writeFrequencies(string_view file_path, const ummss& mapping) {
    cerr << "Writing " << file_path << "\n";
    ofstream report(file_path.data());

    if (!report.is_open()) {
        throw runtime_error("Problem opening frequency file.\n");
    }

    report << "OBJECT\tFREQUENCY\n";
    for (const auto& entry : mapping) {
        report << entry.first << "\t" << mapping.count(entry.first);
    }
}

void writeFrequencies(string_view file_path, const ummss& mapping, const vs& keys) {
    cerr << "Writing " << file_path << "\n";
    ofstream report(file_path.data());
    uss unique_keys(keys.begin(), keys.end());

    if (!report.is_open()) {
        throw runtime_error("Problem opening frequency file.\n");
    }

    report << "OBJECT\tFREQUENCY\n";
    for (const auto& value : keys) {
        report << value << "\t" << mapping.count(value);
    }
}


void writeMeasures(ofstream& report, const measures_result& measures, string_view label1, string_view label2) {
    if (!report.is_open()) {
        throw runtime_error("Could not write measures to report.");
    }
    report << "Average " << label1 << " per " << label2 << ": " << measures.avg << "\n";
    report << "Min: " << measures.min << "\n";
    report << "Max: " << measures.max << "\n";
}

void writeMeasures(ofstream& report, const ummss& mapping, string_view label1, string_view label2) {
    if (!report.is_open()) {
        throw runtime_error("Could not write measures to report.");
    }
    auto measures = calculateMeasures(mapping);
    writeMeasures(report, measures, label1, label2);
}

double getJaccardSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2) {
    if(set1.count() == 0 && set2.count() == 0){
        return 1.0;
    } else{
        double intersection_size = (set1 & set2).count();
        double union_size = (set1 | set2).count();
        return intersection_size / union_size;
    }
}

double getOverlapSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2) {
    if(set1.count() == 0 || set2.count() == 0){
        return 1.0;
    } else{
        double intersection_size = (set1 & set2).count();
        return intersection_size / min(set1.count(), set2.count());
    }
}

std::vector<double> getScores(const msb &sets, std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>)> score) {

    std::vector<double> result;
    for (auto it1 = sets.begin(); it1 != sets.end(); it1++){
        for (auto it2 = next(it1, 1); it2 != sets.end(); it2++){
            result.push_back(score(it1->second, it2->second));
//            std::cerr << "Score (" << it1->first << ", " << it2->first << ") : " << *result.rbegin() << "\n";
        }
    }

    return result;
}

void writeScores(const msb &sets, std::vector<double> scores, std::string_view path_output) {
    ofstream output(path_output.data());
    if (!output.is_open()) {
        throw runtime_error("Problem opening scores file.\n");
    }
    int score = 0;
    output << "SET1\tSET2\tSCORE\n";
    for (auto it1 = sets.begin(); it1 != sets.end(); it1++){
        for (auto it2 = next(it1, 1); it2 != sets.end(); it2++){
            output << it1->first << "\t" << it2->first << "\t" << scores[score] << "\n";
            score++;
        }
    }
    output.close();
}
