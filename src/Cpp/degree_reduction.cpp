#include "degree_reduction.hpp"

using namespace std;

namespace degree_reduction {

/* Requirements:
-- The dataset must contain the mapping for genes, proteins and proteoforms to reactions and pathways.
-- The gene mapping file should have the mapping from genes to proteins.
-- The protein mapping file should have the pathway and reaction names.*/
void doAnalysis(pathway::dataset pathwayDataSet, string_view report_file_path) {
    // Calculate average number of reactions and pathways where a gene, protein and proteoform participate
    double avg_reactions_genes;
    double avg_pathways_genes;
    double avg_reactions_proteins;
    double avg_pathways_proteins;
    double avg_reactions_proteoforms;
    double avg_pathways_proteoforms;

    

    // Calculate average number of proteoforms for a protein
}
   
}  // namespace degree_reduction
