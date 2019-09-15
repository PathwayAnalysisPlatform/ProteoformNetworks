#ifndef MODIFIED_OVERLAP_H_
#define MODIFIED_OVERLAP_H_

#include <charconv>
#include <regex>
#include <fstream>

#include "../overlap_analysis.hpp"
#include "proteoform.hpp"
#include "reactome.hpp"

const double MIN_MODIFIED_ALL_MEMBERS_RATIO = 0.1;
const double MIN_MODIFIED_OVERLAP_MEMBERS_RATIO = 0.9;

namespace modified_overlap {
	void doAnalysis(std::string_view path_file_gene_search,
		std::string_view path_file_protein_search,
		std::string_view path_file_proteoform_search,
		std::string_view path_file_PheGenI,
		std::string_view path_file_mapping_proteins_to_genes,
		std::string_view path_file_report_pathway,
		std::string_view path_file_modified_overlap_pathway_proteins,
		std::string_view path_file_modified_overlap_pathway_proteoforms,
		std::string_view path_file_modified_overlap_pathway_modifications,
		std::string_view path_file_report_trait,
		std::string_view path_file_modified_overlap_trait_proteins,
		std::string_view path_file_modified_overlap_trait_proteoforms,
		std::string_view path_file_modified_overlap_trait_modifications);

	template <size_t total_proteoforms>
	void writeReportRecords(
		std::ofstream& output,
		const std::map<std::pair<std::string, std::string>, std::bitset<total_proteoforms>>& examples,
		const umss& pathways_to_names,
		const reactome_proteoform_sets& proteoform_pathways,
		const bimap_str_int& proteoforms) {
		for (const auto& example : examples) {
			output << example.first.first << "\t" << example.first.second << "\t" << pathways_to_names.at(example.first.first) << "\t"
				<< pathways_to_names.at(example.first.second) << "\t";
			output << proteoform_pathways.at(example.first.first).count() << "\t" << proteoform_pathways.at(example.first.second).count() << "\t";
			output << example.second.count() << "\t";
			printMembers(output, example.second, proteoforms);
			output << "\n";
		}
	}
}

#endif /* MODIFIED_OVERLAP_H_ */