#include "percolation.hpp"

using namespace std;

namespace percolation {

void doAnalysis(const string& path_file_gene_search,
                const string& path_file_protein_search,
                const string& path_file_proteoform_search) {
   const auto [adjacency_list_proteins, adjacency_list_proteoforms] = loadReactomeNetworks(path_file_gene_search, path_file_protein_search, path_file_proteoform_search);

    

   // TODO: Sample to create approximated link and node percolation curve for pp, pm and mm
   // TODO: Plot and save files of percolation curve.
   // TODO: Calculate the theoretical percolation threshold and compare it with the approximated curves.

   // TODO: Extract subgroups and create new percolation curves.

   // TODO: Get most common PTM types for each full network and the subgroups.
   // TODO: Get most common Pathways for the members of each full network and the subgroups.
}

}  // namespace percolation