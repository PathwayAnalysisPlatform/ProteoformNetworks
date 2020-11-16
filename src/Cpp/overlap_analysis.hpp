#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <iostream>
#include <ctime>
#include <utility>

#include "scores.hpp"
#include "bimap_str_int.hpp"
#include "others/uniprot.hpp"

#include "Interactome.hpp"
#include "create_modules.hpp"

#include <windows.h>
#include "tools/files.hpp"

int countSharedVertices(Module m1, Module m2, Interactome interactome);

double getOverlapSize(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

double getOverlapSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

// Calculate Jaccard index, which is intersection over union
double getJaccardSimilarity(base::dynamic_bitset<> set1, base::dynamic_bitset<> set2);

void calculateOverlap(Interactome interactome, std::vector<std::map<std::string, Module>> modules,
                      const std::string &output_path);

#endif /* OVERLAP_H_ */