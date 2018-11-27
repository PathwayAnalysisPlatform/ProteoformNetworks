#ifndef ARTEFACTUAL_OVERLAP_H_
#define ARTEFACTUAL_OVERLAP_H_

#include "../overlap.hpp"

namespace artefactual_overlap {
void doAnalysis(const std::string& path_file_gene_search,
                const std::string& path_file_protein_search,
                const std::string& path_file_proteoform_search,
                const std::string& path_file_PheGenI_full,
                const std::string& path_file_mapping_proteins_to_genes,
                const std::string& path_file_report_pathway,
                const std::string& path_file_report_trait);
}

#endif /* ARTEFACTUAL_OVERLAP_H_ */