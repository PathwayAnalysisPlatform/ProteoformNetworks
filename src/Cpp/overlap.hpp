#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <iostream>

#include "phegeni.hpp"
#include "bimap.hpp"

const int MIN_OVERLAP_SIZE = 1;
const int MAX_OVERLAP_SIZE = 100;

const int MIN_SET_SIZE = 1;
const int MAX_SET_SIZE = 200;

void doOverlapAnalysis(std::string_view path_file_PheGenI, std::string_view path_file_gene_search);

#endif /* OVERLAP_H_ */