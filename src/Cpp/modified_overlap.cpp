#include "modified_overlap.hpp"

using namespace std;

namespace modified_overlap {

	void writePathwayReport(std::ofstream& output,
		const std::map<std::pair<std::string, std::string>, std::bitset<REACTOME_PROTEOFORMS>>& examples,
		const umss& pathways_to_names,
		const um < std::string, std::bitset<REACTOME_PROTEOFORMS >>& pathways_to_proteoforms,
		const bimap& proteoforms) {
		output << "PATHWAY_1\tPATHWAY_2\tPATHWAY_1_NAME\tPATHWAY_2_NAME\t";
		output << "PATHWAY_1_PROTEOFORM_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\tOVERLAP_SIZE\t";
		output << "OVERLAP_PROTEOFORMS\n";
		for (const auto& example : examples) {
			output << example.first.first << "\t" << example.first.second << "\t" << pathways_to_names.at(example.first.first) << "\t"
				<< pathways_to_names.at(example.first.second) << "\t";
			output << pathways_to_proteoforms.at(example.first.first).count() << "\t" << pathways_to_proteoforms.at(example.first.second).count() << "\t";
			output << example.second.count() << "\t";
			printMembers(output, example.second, proteoforms);
			output << "\n";
		}
	}

	void writePhenotypeReport(std::ofstream& output,
		const map<pair<string, string>, bitset<PHEGENI_PROTEOFORMS>>& examples,
		const umss& pathways_to_names,
		const um<string, bitset<PHEGENI_PROTEOFORMS>>& pathways_to_proteoforms,
		const bimap& proteoforms) {
		output << "PHENOTYPE_1\tPHENOTYPE_1_2\tPHENOTYPE_1_NAME\tPHENOTYPE_2_NAME\t";
		output << "PHENOTYPE_1_PROTEOFORM_SIZE\tPHENOTYPE_2_PROTEOFORM_SIZE\tOVERLAP_SIZE\t";
		output << "OVERLAP_PROTEOFORMS\n";
		for (const auto& example : examples) {
			output << example.first.first << "\t" << example.first.second << "\t" << pathways_to_names.at(example.first.first) << "\t"
				<< pathways_to_names.at(example.first.second) << "\t";
			output << pathways_to_proteoforms.at(example.first.first).count() << "\t" << pathways_to_proteoforms.at(example.first.second).count() << "\t";
			output << example.second.count() << "\t";
			printMembers(output, example.second, proteoforms);
			output << "\n";
		}
	}

	void writeFrequencies(std::string_view modifications_file_path, std::string_view proteins_file_path, std::string_view proteoforms_file_path,
		const map<string, int>& mod_to_freq, const map<string, int>& proteins_to_freq, const map<string, int>& proteoforms_to_freq) {
		std::ofstream modifications_file(modifications_file_path.data());
		std::ofstream proteins_file(proteins_file_path.data());
		std::ofstream proteoforms_file(proteoforms_file_path.data());

		modifications_file << "MODIFICATION\tFREQUENCY\n";
		for (const auto& modification : mod_to_freq) {
			modifications_file << modification.first << "\t" << modification.second << "\n";
		}
		cerr << "Finished writing modification frequencies.\n";

		proteins_file << "PROTEIN\tFREQUENCY\n";
		for (const auto& protein : proteins_to_freq) {
			proteins_file << protein.first << "\t" << protein.second << "\n";
		}
		cerr << "Finished writing protein frequencies.\n";

		proteoforms_file << "PROTEOFORM\tFREQUENCY\n";
		for (const auto& proteoform : proteoforms_to_freq) {
			proteoforms_file << proteoform.first << "\t" << proteoform.second << "\n";
		}
		cerr << "Finished writing proteoform frequencies.\n";
	}

	template <size_t total_proteoforms>
	Frequencies getFrequencies(const map<pair<string, string>, bitset<total_proteoforms>>& proteoform_overlap_pairs,
		const bimap& proteoforms) {
		Frequencies frequencies;
		for (const auto& overlap_pair : proteoform_overlap_pairs) {
			// For each overlap set
			for (int I = 0; I < overlap_pair.second.size(); I++) {
				if (overlap_pair.second.test(I)) {
					string accession = proteoform::getAccession(proteoforms.entities[I]);
					if (frequencies.proteins.find(accession) == frequencies.proteins.end()) {
						frequencies.proteins.emplace(accession, 0);
					}
					frequencies.proteins[accession]++;

					if (frequencies.proteoforms.find(proteoforms.entities[I]) == frequencies.proteoforms.end()) {
						frequencies.proteoforms.emplace(proteoforms.entities[I], 0);
					}
					frequencies.proteoforms[proteoforms.entities[I]]++;

					for (const auto& modification : proteoform::getModifications(proteoforms.entities[I])) {
						if (frequencies.modifications.find(modification) == frequencies.modifications.end()) {
							frequencies.modifications[modification] = 0;
						}
						frequencies.modifications[modification]++;
					}
				}
			}
		}
		return frequencies;
	}

	void plotFrequencies(const std::string& report_file_path, const std::string& modifications_file_path, const std::string& proteins_file_path, const std::string& proteoforms_file_path) {
		string command = "Rscript src/R/modified_overlap.R " + report_file_path + " " +
			modifications_file_path + " " + proteins_file_path + " " + proteoforms_file_path;
		cerr << "Plotting modified overlap frequencies: " << command << endl;
		system(command.c_str());
	}

	void reportPathwayPairs(std::string_view path_file_proteoform_search,
		std::string_view report_file_path,
		std::string_view modifications_file_path,
		std::string_view proteins_file_path,
		std::string_view proteoforms_file_path) {
		const auto proteoforms = createBimap(path_file_proteoform_search);
		const auto pathways_to_names = loadPathwayNames(path_file_proteoform_search);
		const auto pathways_to_proteoforms = loadProteoformSets(path_file_proteoform_search, proteoforms, true);
		const bitset<REACTOME_PROTEOFORMS> modified_proteoforms = proteoform::getSetOfModifiedProteoforms<REACTOME_PROTEOFORMS>(proteoforms);

		cout << "Reporting pathway pairs with only modified overlap...\n";

		// Compare all the pairs of selected pathways and select pairs that overlap only in a percentage of modified proteins
		const auto& examples = findOverlappingProteoformSets(pathways_to_proteoforms, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE,
			modified_proteoforms, MIN_MODIFIED_ALL_MEMBERS_RATIO, MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

		std::ofstream report(report_file_path.data());
		writePathwayReport(report, examples, pathways_to_names, pathways_to_proteoforms, proteoforms);

		const auto [mod_to_freq, proteins_to_freq, proteoforms_to_freq] = getFrequencies(examples, proteoforms);
		writeFrequencies(modifications_file_path, proteins_file_path, proteoforms_file_path, mod_to_freq, proteins_to_freq, proteoforms_to_freq);
		plotFrequencies(report_file_path.data(), modifications_file_path.data(), proteins_file_path.data(), proteoforms_file_path.data());
	}

	void reportPhenotypePairs(std::string_view path_file_gene_search,
		std::string_view path_file_protein_search,
		std::string_view path_file_proteoform_search,
		std::string_view path_file_PheGenI,
		std::string_view path_file_mapping_proteins_genes,
		std::string_view report_file_path,
		std::string_view modifications_file_path,
		std::string_view proteins_file_path,
		std::string_view proteoforms_file_path) {
		// Load data: trait gene sets, trait proteoform sets
		cout << "Loading PheGen data\n";
		const auto reactome_genes = createBimap(path_file_gene_search);
		const auto [phegeni_genes, phegeni_traits] = loadPheGenIEntities(path_file_PheGenI, reactome_genes);
		const auto [genes_to_proteins, proteins_to_genes] = loadMappingGenesProteins(path_file_mapping_proteins_genes.data());
		const auto [proteins_to_proteoforms, proteoforms_to_proteins] = loadMappingProteinsProteoforms(path_file_proteoform_search.data());
		const auto phegeni_proteins = deductProteinsFromGenes(path_file_mapping_proteins_genes, genes_to_proteins, phegeni_genes);
		const auto phegeni_proteoforms = deductProteoformsFromProteins(proteins_to_proteoforms, phegeni_proteins);
		const bitset<PHEGENI_PROTEOFORMS> modified_proteoforms = proteoform::getSetOfModifiedProteoforms<PHEGENI_PROTEOFORMS>(phegeni_proteoforms);

		const auto [adjacency_list_proteins, adjacency_list_proteoforms] = loadReactomeNetworks(path_file_gene_search, path_file_protein_search, path_file_proteoform_search);
		const auto [traits_to_genes, genes_to_traits] = loadPheGenISets(path_file_PheGenI.data(), reactome_genes, phegeni_genes, phegeni_traits);
		const auto sets_to_names = createTraitNames(traits_to_genes);
		const auto traits_to_proteins = convertGeneSets(traits_to_genes, phegeni_genes, genes_to_proteins, phegeni_proteins, adjacency_list_proteins);
		const auto traits_to_proteoforms = convertProteinSets(traits_to_proteins, phegeni_proteins, proteins_to_proteoforms, phegeni_proteoforms, adjacency_list_proteoforms);

		// Calculate overlap
		const auto overlapping_proteoform_set_pairs = findOverlappingPairs(traits_to_proteoforms, 1, PHEGENI_PROTEOFORMS, 1, PHEGENI_PROTEOFORMS);
		const auto non_overlapping_proteoform_set_pairs = findOverlappingPairs(traits_to_proteoforms, 0, 0, 1, PHEGENI_PROTEOFORMS);
		const auto examples = findOverlappingProteoformSets(traits_to_proteoforms, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE,
			modified_proteoforms, MIN_MODIFIED_ALL_MEMBERS_RATIO, MIN_MODIFIED_OVERLAP_MEMBERS_RATIO);

		// Write report
		std::ofstream report(report_file_path.data());
		writePhenotypeReport(report, examples, sets_to_names, traits_to_proteoforms, phegeni_proteoforms);

		// Analyse modifications
		const auto [mod_to_freq, proteins_to_freq, proteoforms_to_freq] = getFrequencies(examples, phegeni_proteoforms);
		writeFrequencies(modifications_file_path, proteins_file_path, proteoforms_file_path, mod_to_freq, proteins_to_freq, proteoforms_to_freq);
		plotFrequencies(report_file_path.data(), modifications_file_path.data(), proteins_file_path.data(), proteoforms_file_path.data());
	}

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
		std::string_view path_file_modified_overlap_trait_modifications,
		std::string_view path_file_modified_overlap_trait_proteins,
		std::string_view path_file_modified_overlap_trait_proteoforms) {
		cout << "Searching for modified overlap examples...\n";

		// Part 1: Find examples of pathways that overlap only in modified proteins
		reportPathwayPairs(path_file_proteoform_search,
			path_file_report_pathway,
			path_file_modified_overlap_pathway_modifications,
			path_file_modified_overlap_pathway_proteins,
			path_file_modified_overlap_pathway_proteoforms);

		// Part 2: Find examples of disease modules that overlap only in modified proteins
		reportPhenotypePairs(path_file_gene_search,
			path_file_protein_search,
			path_file_proteoform_search,
			path_file_PheGenI,
			path_file_mapping_proteins_to_genes,
			path_file_report_trait,
			path_file_modified_overlap_trait_modifications,
			path_file_modified_overlap_trait_proteins,
			path_file_modified_overlap_trait_proteoforms);
	}

}  // namespace modified_overlap
