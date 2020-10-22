#include "overlap_analysis.hpp"

double calculateJaccardIndex(Module m1, Module m2);

double calculateOverlapCoefficient(Module m1, Module m2);

int countSharedVertices(Module m1, Module m2);

using namespace std;

std::string GetExeFileName() {
    char buffer[MAX_PATH];
    GetModuleFileName(NULL, buffer, MAX_PATH);
    return std::string(buffer);
}

std::string GetExePath() {
    std::string f = GetExeFileName();
    return f.substr(0, f.find_last_of("\\/"));
}


int countSharedVertices(Module m1, Module m2, Interactome interactome) {

    int count = 0;

    for(auto& entry : m1.adj){
        int v = entry.first;
        if(!interactome.isSimpleEntity(v)){
           if(m2.hasVertex(v))
               count++;
        }
    }

    return count;
}

double getOverlapSize(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2) {
    return (set1 & set2).count();
}

double getOverlapSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2) {
    if (set1.count() == 0 || set2.count() == 0) {
        return 1.0;
    } else {
        double intersection_size = (set1 & set2).count();
        return intersection_size / min(set1.count(), set2.count());
    }
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

// Create 3 column table: Module1, Module2, Level, # Shared vertices, Overlap Coefficient, Jaccard index
void calculateOverlap(Interactome interactome, std::vector<std::map<std::string, Module>> modules,
                      const std::string &output_path) {

    std::ofstream fall(output_path + "all_pairs.tsv");
    std::ofstream foverlapping(output_path + "overlapping_pairs.tsv");

    if (!fall.is_open()) {
        std::string message = "Cannot open overlap results file " + output_path + "all_pairs.tsv";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    if (!foverlapping.is_open()) {
        std::string message = "Cannot open overlap results file " + output_path + "overlapping_pairs.tsv";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    fall << "LEVEL\t" << "MODULE1\t" << "MODULE2\t" << "SHARED_ACCESSIONED_ENTITIES\t" << "OVERLAP_COEFFICIENT\t" << "JACCARD_INDEX\n";
    foverlapping << "LEVEL\t" << "MODULE1\t" << "MODULE2\t" << "SHARED_ACCESSIONED_ENTITIES\t" << "OVERLAP_COEFFICIENT\t" << "JACCARD_INDEX\n";

    for (int level = 0; level < 3; level++) {
        std::cerr << "Calculating overlap scores for " << LEVELS[level] << std::endl;

        int cont = 1;
        for(auto it1 = modules[level].begin(); it1 != modules[level].end(); it1++){
            cout << cont << ": " << it1->first << "\n";
            auto it2 = it1;
            for(it2++; it2 != modules[level].end(); it2++) {
//                cout << " -- with " << it2->first << std::endl;
                int sharedVertices = getOverlapSize(it1->second.accessioned_entity_vertices, it2->second.accessioned_entity_vertices);
                double overlapCoefficient = 0;
                double jaccardIndex = 0;

                fall << LEVELS[level] << "\t" << it1->first << "\t" << it2->first << "\t";
                fall << sharedVertices << "\t" << overlapCoefficient << "\t" << jaccardIndex << "\n";

                if(sharedVertices){

                    double overlapCoefficient = getOverlapSimilarity(it1->second.accessioned_entity_vertices, it2->second.accessioned_entity_vertices);
                    double jaccardIndex = getJaccardSimilarity(it1->second.accessioned_entity_vertices, it2->second.accessioned_entity_vertices);

                    foverlapping << LEVELS[level] << "\t" << it1->first << "\t" << it2->first << "\t";
                    foverlapping << sharedVertices << "\t" << overlapCoefficient << "\t" << jaccardIndex << std::endl;
                }
            }
            cont++;
        }
    }
}


int main(int argc, char *argv[]) try {

    if (argc < 4) {
        std::cerr << "Missing arguments. Expected: 8 arguments:\n\n"
                  << " * - [1] File with Gene sets to create disease modules\n"
                  << " * - [2] File of Interactome vertices\n"
                  << " * - [3] File of Interactome edges (with indexed vertices)\n"
                  << " * - [4] File with ranges: start and end vertex index for each entity type\n"
                  << " * - [5] Mapping from proteins to genes\n"
                  << " * - [6] Mapping from proteins to proteoforms\n"
                  << " * - [7] Modules output path\n"
                  << " * - [8] Overlap output path";
        throw std::runtime_error("Missing arguments.");
        return 0;
    }

    std::string file_phegeni = argv[1];
    std::string file_vertices = argv[2];
    std::string file_edges = argv[3];
    std::string file_ranges = argv[4];
    std::string file_proteins_to_genes = argv[5];
    std::string file_proteins_to_proteoforms = argv[6];
    std::string modules_output_path = argv[7];
    std::string overlap_output_path = argv[8];

    std::cout << "Reading Interactome...\n\n";
    Interactome interactome(file_vertices, file_edges, file_ranges, file_proteins_to_genes,
                            file_proteins_to_proteoforms);

    std::cout << "Creating disease modules.\n\n";
    auto modules = createModules(file_phegeni, interactome, modules_output_path);

    std::cout << "Calculate overlap.\n\n";
    calculateOverlap(interactome, modules, overlap_output_path);
}
catch (const std::exception &ex) {
    std::cout << ex.what() << "\n";
}
