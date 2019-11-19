#include "scores.hpp"
#include <iomanip>

using namespace std;

const measures_result calculateMeasures(const ummss &mapping) {
    double min = static_cast<double>(mapping.count(mapping.begin()->first));
    double max = static_cast<double>(mapping.count(mapping.begin()->first));
    double avg = 0.0;
    double sum = 0.0;
    for (const auto &entry : mapping) {
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
const measures_result calculateMeasuresWithSelectedKeys(const ummss &mapping,
                                                        const vs &keys) {
    // Calculate averag proteoforms for all them
    double min = static_cast<double>(mapping.count(*keys.begin()));
    double max = static_cast<double>(mapping.count(*keys.begin()));
    double avg = 0.0;
    double sum = 0.0;
    for (const auto &value : keys) {
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

void writeFrequencies(ofstream &report, const ummss &mapping, string_view label) {
    if (!report.is_open()) {
        throw runtime_error("Problem opening frequency file.\n");
    }

    for (const auto &entry : mapping) {
        report << entry.first << "\t" << mapping.count(entry.first);
    }
}

void writeFrequencies(string_view file_path, const ummss &mapping) {
    cerr << "Writing " << file_path << "\n";
    ofstream report(file_path.data());

    if (!report.is_open()) {
        throw runtime_error("Problem opening frequency file.\n");
    }

    report << "OBJECT\tFREQUENCY\n";
    for (const auto &entry : mapping) {
        report << entry.first << "\t" << mapping.count(entry.first);
    }
}

void writeFrequencies(string_view file_path, const ummss &mapping, const vs &keys) {
    cerr << "Writing " << file_path << "\n";
    ofstream report(file_path.data());
    uss unique_keys(keys.begin(), keys.end());

    if (!report.is_open()) {
        throw runtime_error("Problem opening frequency file.\n");
    }

    report << "OBJECT\tFREQUENCY\n";
    for (const auto &value : keys) {
        report << value << "\t" << mapping.count(value);
    }
}


void writeMeasures(ofstream &report, const measures_result &measures, string_view label1, string_view label2) {
    if (!report.is_open()) {
        throw runtime_error("Could not write measures to report.");
    }
    report << "Average " << label1 << " per " << label2 << ": " << measures.avg << "\n";
    report << "Min: " << measures.min << "\n";
    report << "Max: " << measures.max << "\n";
}

void writeMeasures(ofstream &report, const ummss &mapping, string_view label1, string_view label2) {
    if (!report.is_open()) {
        throw runtime_error("Could not write measures to report.");
    }
    auto measures = calculateMeasures(mapping);
    writeMeasures(report, measures, label1, label2);
}

double getJaccardSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2) {
    if (set1.count() == 0 && set2.count() == 0) {
        return 1.0;
    } else {
        double intersection_size = (set1 & set2).count();
        double union_size = (set1 | set2).count();
        return intersection_size / union_size;
    }
}

double getOverlapSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2) {
    if (set1.count() == 0 || set2.count() == 0) {
        return 1.0;
    } else {
        double intersection_size = (set1 & set2).count();
        return intersection_size / min(set1.count(), set2.count());
    }
}

double getOverlapSize(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2) {
    return (set1 & set2).count();
}

pair_map<double> getScores(const vb &sets,
                           std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>)> score_function) {

    pair_map<double> result;
    for (int I1 = 0; I1 < sets.size(); I1++) {
        for (int I2 = I1 + 1; I2 < sets.size(); I2++) {
            result[make_pair(I1, I2)] = score_function(sets[I1], sets[I2]);
//            std::cerr << "Score (" << it1->first << ", " << it2->first << ") : " << *result.rbegin() << "\n";
        }
    }

    return result;
}

void writeScores(const bimap_str_int &groups,
                 const modules &entity_modules,
                 const std::vector<std::string> &features_labels,
                 const std::vector<pair_map<double>> &features,
                 std::string_view file_output) {
    ofstream output(file_output.data());
    if (!output.is_open()) {
        throw runtime_error("Problem opening scores file.\n");
    }
    int score = 0;
    output << "TRAIT1\tTRAIT2";
    for (const auto &label : features_labels)
        output << "\t" << label;
    output << "\n";
    for (int I1 = 0; I1 < entity_modules.group_to_members.size(); I1++) {
        for (int I2 = I1 + 1; I2 < entity_modules.group_to_members.size(); I2++) {
            const std::pair<int, int> index_pair = make_pair(I1, I2);
            output << groups.int_to_str[I1] << "\t" << groups.int_to_str[I2] << std::setprecision(5);
            for (const auto &feature : features)
                output << "\t" << feature.at(index_pair);
            output << "\n";
            score++;
        }
    }
    output.close();
}


