#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <iostream>
#include <ctime>
#include <sys/stat.h>
#include <utility>

#include "scores.hpp"
#include "bimap_str_int.hpp"
#include "others/uniprot.hpp"

#include "Interactome.hpp"
#include "create_modules.hpp"

#include <windows.h>

inline bool file_exists(const std::string &name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

double getOverlapSize(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

double getOverlapSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

// Calculate Jaccard index, which is intersection over union
double getJaccardSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

void calculateOverlap(Interactome interactome, std::vector<std::map<std::string, Module>> modules, const std::string& output_path);


#endif /* OVERLAP_H_ */