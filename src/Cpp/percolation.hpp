#ifndef PERCOLATION_H_
#define PERCOLATION_H_

#include "overlap.hpp"

namespace percolation {

void doAnalysis(std::string_view  path_file_gene_search,
                std::string_view  path_file_protein_search,
                std::string_view  path_file_proteoform_search);
}
#endif /* PERCOLATION_H_ */