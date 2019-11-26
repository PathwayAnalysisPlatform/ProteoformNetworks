#include "scores.hpp"


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

/*
 * Calculates the score for each pair of sets.
 * The calculation requires just the sets of vertices. No need for the edges.
 */
pair_map<double> getScores(const vb &vertex_sets,
                           std::function<double(base::dynamic_bitset<>,
                                                base::dynamic_bitset<>)> score_function,
                           const int min_module_size, const int max_module_size) {

    pair_map<double> result;
    for (int I1 = 0; I1 < vertex_sets.size(); I1++) {
        if (vertex_sets[I1].count() < min_module_size || vertex_sets[I1].count() > max_module_size) continue;
        for (int I2 = I1 + 1; I2 < vertex_sets.size(); I2++) {
            if (vertex_sets[I2].count() < min_module_size || vertex_sets[I2].count() > max_module_size) continue;
            std::cerr << I1 << " -- " << I2 << std::endl;
            auto score = score_function(vertex_sets[I1], vertex_sets[I2]);
            if (score > 0)
                result[make_pair(I1, I2)] = score;
        }
    }

    return result;
}

pair_map<double> getScores(const vb &vertex_sets,
                           std::function<double(base::dynamic_bitset<>,
                                                base::dynamic_bitset<>)> score_function,
                           const pair_map<double> &prev_scores) {

    pair_map<double> result;
    for (const auto &pair : prev_scores) {
//        std::cerr << "Try " << pair.first.first << " , " << pair.first.second << std::endl;
        result[pair.first] = score_function(vertex_sets[pair.first.first], vertex_sets[pair.first.second]);
//        std::cerr << "( " << pair.first.first << " , " << pair.first.second << " ) = " << result[pair.first]
//                  << std::endl;
    }

    std::cerr << "Finished calculating scores. " << std::endl;
    return result;
}

/*
 * Calculates the score for each pair of sets.
 * The calculation requires the vertices and edges represented by the sets.
 */
pair_map<double>
getScores(const vb &vertex_sets,
          const vusi &edges,
          std::function<double(base::dynamic_bitset<>, base::dynamic_bitset<>, vusi)> score_function,
          const pair_map<double> &overlap_sizes) {

    pair_map<double> result;
    for (const auto &pair : overlap_sizes) {
        result[pair.first] = score_function(vertex_sets[pair.first.first], vertex_sets[pair.first.second], edges);
        std::cerr << "( " << pair.first.first << " , " << pair.first.second << " ) = " << result[pair.first]
                  << std::endl;
    }

    return result;
}

/*
 * Calculates the number of nodes in the interface between the two modules
 */
double calculate_interface_size_nodes(const base::dynamic_bitset<> &V1,
                                      const base::dynamic_bitset<> &V2,
                                      const vusi &E) {
    std::set<int> interface_nodes;
    // Count the number of nodes with an edge going to the other set
    for (int start_node = 0; start_node < V1.size(); start_node++) {
        for (int end_node : E[start_node]) {
            if (V1[start_node] && V2[end_node]) {
                interface_nodes.insert({start_node, end_node});
            } else if (V1[start_node] && V2[start_node]) {
                interface_nodes.insert(start_node);
            } else if (V1[end_node] && V2[end_node]) {
                interface_nodes.insert(end_node);
            }
        }
    }
    return interface_nodes.size();
}

double calculate_interface_size_edges(const base::dynamic_bitset<> &V1,
                                      const base::dynamic_bitset<> &V2,
                                      const vusi &E) {
    int result = 0;

    for (int start_node = 0; start_node < V1.size(); start_node++) {
        if (V1[start_node] || V2[start_node]) {             // Check if the node belongs to any of the two modules
            for (int end_node : E[start_node]) {
                if (start_node < end_node) {                // Check all edges only once, not in both directions.
                    if (V1[start_node] || V2[start_node]) { // Check if the node belongs to any of the two modules
                        // Check if the nodes are both in the overlap Or if the nodes are in different modules
                        // Or one node is in the overlap and the other not
                        if (!((!V1[start_node] && !V1[end_node] && V2[start_node] && V2[end_node]) ||
                              (V1[start_node] && V1[end_node] && !V2[start_node] && !V2[end_node]))) {
                            result++;
                            std::cerr << "( " << start_node << ", " << end_node << ")\n";
                        }
                    }
                }
            }
        }
    }
    return result;
}

pair_map<double> getScores(const vb &sets, std::function<double()> score_function);

void writeScores(const bimap_str_int &groups,
                 const modules &entity_modules,
                 const std::vector<std::string> &features_labels,
                 const std::vector<pair_map<double>> &features,
                 std::string_view file_output) {
    assert(features.size() > 0);

    ofstream output(file_output.data());
    if (!output.is_open()) {
        throw runtime_error("Problem opening scores file.\n");
    }
    int score = 0;
    output << "TRAIT1\tTRAIT2";
    for (const auto &label : features_labels)
        output << "\t" << label;
    output << "\n";
    for (const auto &pair : features[0]) {
        output << groups.int_to_str[pair.first.first] << "\t" << groups.int_to_str[pair.first.second] << std::setprecision(5);
        for (const auto &feature : features)
            output << "\t" << feature.at(pair.first);
        output << "\n";
        score++;
    }
    output.close();
}


