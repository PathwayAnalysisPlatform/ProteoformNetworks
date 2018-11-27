#ifndef PERCOLATION_H_
#define PERCOLATION_H_

#include "../overlap.hpp"

namespace percolation {

void doAnalysis(const std::string& path_file_gene_search,
                const std::string& path_file_protein_search,
                const std::string& path_file_proteoform_search);
}
#endif /* PERCOLATION_H_ */