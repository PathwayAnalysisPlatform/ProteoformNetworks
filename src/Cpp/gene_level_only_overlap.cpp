#include "gene_level_only_overlap.hpp"

using namespace std;

namespace gene_level_only_overlap {

	/**
	 * Find gene level only overlaps: pairs of pathways that share nodes only in the
	 * gene or protein level, but not at the proteoform level
	 */
	template <size_t size_genes, size_t size_proteins, size_t size_proteoforms>
	set<pair<string, string>> findGeneLevelOnlyOverlapPairs(const um<string, bitset<size_proteoforms>>& sets_to_proteoforms,
		const map<pair<string, string>, bitset<size_genes>>& overlapping_gene_set_pairs,
		const map<pair<string, string>, bitset<size_proteins>>& overlapping_protein_set_pairs,
		const map<pair<string, string>, bitset<size_proteoforms>>& non_overlapping_proteoform_set_pairs) {
		set<pair<string, string>> result;

		cerr << "Comparing gene and proteoform pairs..." << endl;
		for (const auto& gene_set_pair : overlapping_gene_set_pairs) {                                                          // For each pair of overlapping sets at gene level
			if (non_overlapping_proteoform_set_pairs.find(gene_set_pair.first) != non_overlapping_proteoform_set_pairs.end()) {  // Check if they do not overlap
				result.insert(gene_set_pair.first);
			}
		}
		cerr << "Comparing protein and proteoform pairs..." << endl;
		for (const auto& protein_set_pair : overlapping_protein_set_pairs) {
			if (non_overlapping_proteoform_set_pairs.find(protein_set_pair.first) != non_overlapping_proteoform_set_pairs.end()) {
				result.insert(protein_set_pair.first);
			}
		}
		cout << "Found " << result.size() << " examples.\n";
		return result;
	}

	// Select the proteoforms in the set which have any of the argument accessions
	template <size_t size_proteoforms>
	bitset<size_proteoforms> getProteoformSetByAccessionStrings(const uss& accessions,
		const bitset<size_proteoforms>& proteoform_set,
		const bimap& proteoforms) {
		bitset<size_proteoforms> result;
		for (int I = 0; I < proteoform_set.size(); I++) {
			if (proteoform_set.test(I)) {
				if (accessions.find(proteoform::getAccession(proteoforms.entities.at(I))) != accessions.end()) {
					result.set(I);
				}
			}
		}
		return result;
	}

	template <size_t size_genes, size_t size_proteins, size_t size_proteoforms>
	void writeReportRecords(ofstream& output,
		const set<pair<string, string>> examples,
		const umss& sets_to_names,
		const um<string, bitset<size_genes>>& sets_genes,
		const um<string, bitset<size_proteins>>& sets_proteins,
		const um<string, bitset<size_proteoforms>>& sets_proteoforms,
		const bimap& phegeni_genes,
		const bimap& proteins,
		const bimap& proteoforms) {
		cout << "Writing records...\n";

		auto t0 = clock();

		for (const auto& example : examples) {
			bitset<size_genes> overlap_genes = sets_genes.at(example.first) & sets_genes.at(example.second);
			bitset<size_proteins> overlap_proteins = sets_proteins.at(example.first) & sets_proteins.at(example.second);
			bitset<size_proteoforms> overlap_proteoforms = sets_proteoforms.at(example.first) & sets_proteoforms.at(example.second);
			auto protein_overlap_members = getStringSetFromEntityBitset(overlap_proteins, proteins);

			auto decomposed_overlap_proteoforms_1 = getProteoformSetByAccessionStrings(protein_overlap_members, sets_proteoforms.at(example.first), proteoforms);
			auto decomposed_overlap_proteoforms_2 = getProteoformSetByAccessionStrings(protein_overlap_members, sets_proteoforms.at(example.second), proteoforms);

			output << example.first << "\t" << example.second << "\t" << sets_to_names.at(example.first) << "\t" << sets_to_names.at(example.second) << "\t";
			output << sets_genes.at(example.first).count() << "\t" << sets_proteins.at(example.first).count() << "\t";
			output << sets_proteoforms.at(example.first).count() << "\t";
			output << sets_genes.at(example.second).count() << "\t" << sets_proteins.at(example.second).count() << "\t";
			output << sets_proteoforms.at(example.second).count() << "\t";
			output << overlap_genes.count() << "\t" << overlap_proteins.count() << "\t" << overlap_proteoforms.count() << "\t";
			printMembers(output, overlap_genes, phegeni_genes);
			output << "\t";
			printMembers(output, overlap_proteins, proteins);
			output << "\t";
			printMembers(output, overlap_proteoforms, proteoforms);
			output << "\t";

			// Print the proteoforms to which the genes were converted, and that do not show the overlap anymore
			printMembers(output, decomposed_overlap_proteoforms_1, proteoforms);
			output << "\t";
			printMembers(output, decomposed_overlap_proteoforms_2, proteoforms);
			output << "\n";
		}
		auto t1 = clock();
		cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
	}

	void writePathwayReport(
		ofstream& output,
		const set<pair<string, string>> examples,
		const umss& sets_to_names,
		const reactome_gene_sets& sets_genes,
		const reactome_protein_sets& sets_proteins,
		const reactome_proteoform_sets& sets_proteoforms,
		const bimap& phegeni_genes,
		const bimap& proteins,
		const bimap& proteoforms) {
		cerr << "Reporting pathway pairs with gene level only overlap...\n";
		output << "PATHWAY_1\tPATHWAY_2\tPATHWAY_1_NAME\tPATHWAY_2_NAME\t";
		output << "PATHWAY_1_GENE_SIZE\tPATHWAY_1_PROTEIN_SIZE\tPATHWAY_1_PROTEOFORM_SIZE\t";
		output << "PATHWAY_2_GENE_SIZE\tPATHWAY_2_PROTEIN_SIZE\tPATHWAY_2_PROTEOFORM_SIZE\t";
		output << "GENE_OVERLAP\tPROTEIN_OVERLAP\tPROTEOFORM_OVERLAP\t";
		output << "OVERLAP_GENES\tOVERLAP_PROTEINS\tOVERLAP_PROTEOFORMS\tDECOMPOSED_OVERLAP_PROTEOFORMS_1\tDECOMPOSED_OVERLAP_PROTEOFORMS_2\n";
		writeReportRecords(output, examples, sets_to_names,
			sets_genes, sets_proteins, sets_proteoforms, phegeni_genes, proteins, proteoforms);
	}

	void writePhenotypeReport(ofstream& output,
		const set<pair<string, string>> examples,
		const umss& sets_to_names,
		const um<string, bitset<PHEGENI_GENES>>& sets_to_genes,
		const um<string, bitset<PHEGENI_PROTEINS>>& sets_to_proteins,
		const um<string, bitset<PHEGENI_PROTEOFORMS>>& sets_to_proteoforms,
		const bimap& phegeni_genes,
		const bimap& proteins,
		const bimap& proteoforms) {
		output << "PHENOTYPE_1\tPHENOTYPE_2\tPHENOTYPE_1_NAME\tPHENOTYPE_2_NAME\t";
		output << "PHENOTYPE_1_GENE_SIZE\tPHENOTYPE_1_PROTEIN_SIZE\tPHENOTYPE_1_PROTEOFORM_SIZE\t";
		output << "PHENOTYPE_2_GENE_SIZE\tPHENOTYPE_2_PROTEIN_SIZE\tPHENOTYPE_2_PROTEOFORM_SIZE\t";
		output << "GENE_OVERLAP\tPROTEIN_OVERLAP\tPROTEOFORM_OVERLAP\t";
		output << "OVERLAP_GENES\tOVERLAP_PROTEINS\tOVERLAP_PROTEOFORMS\tDECOMPOSED_OVERLAP_PROTEOFORMS_1\tDECOMPOSED_OVERLAP_PROTEOFORMS_2\n";

		writeReportRecords(output, examples, sets_to_names,
			sets_to_genes, sets_to_proteins, sets_to_proteoforms,
			phegeni_genes, proteins, proteoforms);
	}

	template <size_t total_num_entities>
	void writeSetReport(std::string_view path_file_report,
		const map<string, bitset<total_num_entities>>& sets_to_entities,
		const vs& index_to_entities) {
		ofstream file_report(path_file_report);
		cout << "Writing report: " << path_file_report << "\n";
		if (!file_report.is_open()) {
			throw runtime_error(path_file_report.data());
		}

		for (const auto& set_to_entities : sets_to_entities) {
			file_report << set_to_entities.first << "\t";
			for (int I = 0; I < total_num_entities; I++) {
				if (set_to_entities.second.test(I)) {
					file_report << index_to_entities.at(I) << "\t";
				}
			}
			file_report << "\n";
		}
	}

	void writeEntitiesReport(std::string_view path_file_report, const vs& entities) {
		ofstream file_report(path_file_report.data());

		cerr << "Writing report: " << path_file_report << "\n";
		if (!file_report.is_open()) {
			throw runtime_error(path_file_report.data());
		}

		for (const auto& entity : entities) {
			file_report << entity << "\n";
		}
	}

	map<string, bitset<REACTOME_GENES>> loadReactionsGeneMembers(std::string_view file_path, const map<string, int>& entities_to_index) {
		map<string, bitset<REACTOME_GENES>> result;
		ifstream file_search(file_path.data());
		string field, gene, pathway;
		bitset<REACTOME_GENES> empty_set;

		getline(file_search, field);                // Skip csv header line
		while (getline(file_search, gene, '\t')) {  // Read the members of each pathway // Read GENE
			getline(file_search, field, '\t');       // Read UNIPROT
			getline(file_search, field, '\t');       // Read REACTION_STID
			getline(file_search, field, '\t');       // Read REACTION_DISPLAY_NAME
			getline(file_search, pathway, '\t');     // Read PATHWAY_STID
			getline(file_search, field);             // Read PATHWAY_DISPLAY_NAME

			if (result.find(pathway) == result.end()) {
				result.emplace(pathway, empty_set);
			}
			if (entities_to_index.find(gene) != entities_to_index.end()) {
				result[pathway].set(entities_to_index.find(gene)->second);
			}
		}
		return result;
	}

	ummss loadGenesAdjacencyList(std::string_view search_file_path) {
		ummss adjacenty_list;
		const auto phegeni_genes = createBimap(search_file_path);
		const auto reactions_to_entities = loadGeneSets(search_file_path, phegeni_genes, false);

		for (const auto& reaction_entry : reactions_to_entities) {
			vs members;
			for (int I = 0; I < reaction_entry.second.size(); I++) {
				if (reaction_entry.second.test(I)) {
					members.push_back(phegeni_genes.entities[I]);
				}
			}

			for (const auto& one_member : members) {
				for (const auto& other_member : members) {
					adjacenty_list.emplace(one_member, other_member);
				}
			}
		}

		return adjacenty_list;
	}

	void reportPathwayPairs(std::string_view path_file_gene_search,
		std::string_view path_file_protein_search,
		std::string_view path_file_proteoform_search,
		std::string_view report_file_path) {
		const auto reactome_genes = loadReactomeEntities(path_file_gene_search);
		const auto reactome_proteins = loadReactomeEntities(path_file_protein_search);
		const auto reactome_proteoforms = loadReactomeEntities(path_file_proteoform_search);
		const auto pathways_to_names = loadPathwayNames(path_file_proteoform_search);

		reactome_gene_sets pathways_to_genes = loadGeneSets(path_file_gene_search, reactome_genes, true);
		reactome_protein_sets pathways_to_proteins = loadProteinSets(path_file_protein_search, reactome_proteins, true);
		reactome_proteoform_sets pathways_to_proteoforms = loadProteoformSets(path_file_proteoform_search, reactome_proteoforms, true);

		cout << "Calculating gene sets overlap..." << endl;
		const auto overlapping_gene_set_pairs = findOverlappingPairs(pathways_to_genes, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE);
		cout << "Calculating protein sets overlap..." << endl;
		const auto overlapping_protein_set_pairs = findOverlappingPairs(pathways_to_proteins, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE);
		cout << "Calculating proteoform sets overlap..." << endl;
		const auto overlapping_proteoform_set_pairs = findOverlappingPairs(pathways_to_proteoforms, 1, REACTOME_PROTEOFORMS, 1, REACTOME_PROTEOFORMS);
		const auto non_overlapping_proteoform_set_pairs = findOverlappingPairs(pathways_to_proteoforms, 0, 0, 1, REACTOME_PROTEOFORMS);

		cerr << "Total proteoforms sets: " << pathways_to_proteoforms.size() << "\n";
		cerr << "Parejas: " << (pathways_to_proteoforms.size() * pathways_to_proteoforms.size() - pathways_to_proteoforms.size()) / 2 << "\n";
		cerr << "Overlapping: " << overlapping_proteoform_set_pairs.size() << "\n";
		cerr << "Non-overlapping: " << non_overlapping_proteoform_set_pairs.size() << "\n";
		cerr << "Total pairs defined: " << overlapping_proteoform_set_pairs.size() + non_overlapping_proteoform_set_pairs.size() << "\n";

		cout << "Finding examples of gene level only overlap...\n";
		const auto examples = findGeneLevelOnlyOverlapPairs(pathways_to_proteoforms,
			overlapping_gene_set_pairs,
			overlapping_protein_set_pairs,
			non_overlapping_proteoform_set_pairs);

		ofstream report(report_file_path.data());
		writePathwayReport(report, examples, pathways_to_names,
			pathways_to_genes, pathways_to_proteins, pathways_to_proteoforms,
			reactome_genes, reactome_proteins, reactome_proteoforms);
	}

	void reportPhenotypePairs(
		std::string_view path_file_gene_search,
		std::string_view path_file_protein_search,
		std::string_view path_file_proteoform_search,
		std::string_view path_file_PheGenI,
		std::string_view path_file_mapping_proteins_to_genes,
		std::string_view path_file_report_trait) {
		// Load reference network

		// Load phegen sets

		// Convert to protein sets: reference gene network, phegen gene sets, mapping from genes to proteins
		// For each protein, decide if it is connected to the other proteins in the set or not in the reference network
		// == For each gene in the phegen set:
		// == Use the mapping from genes to proteins and find the corresponding proteins.
		// == For each correspoding protein:
		// == Add to the new phegen protein set if at least one of its neighbours in the reference network is in the phegen set
		// Convert to proteoform sets: reference protein network, phegen protein sets, mapping proteins to proteoforms
		// Find overlapping pairs
		// Write report

		// Create gene sets for each trait
		const auto reactome_genes = loadReactomeEntities(path_file_gene_search);
		const auto [phegeni_genes, phegeni_traits] = loadPheGenIGenesAndTraits(path_file_PheGenI, reactome_genes);
		const auto mapping_traits_genes = loadPheGenISets(path_file_PheGenI, reactome_genes, phegeni_genes, phegeni_traits);
		const auto traits_to_names = createTraitNames(mapping_traits_genes.traits_to_genes);

		// Create proteoform sets for each trait
		const auto [adjacency_list_proteins, adjacency_list_proteoforms] = loadReactomeNetworks(path_file_protein_search, path_file_proteoform_search);		

		const auto [genes_to_proteins, proteins_to_genes] = loadMappingGenesProteins(path_file_mapping_proteins_to_genes);
		const auto phegeni_proteins = deductProteinsFromGenes(genes_to_proteins, phegeni_genes);
		const um<string, bitset<PHEGENI_PROTEINS>> traits_to_proteins = convertGeneSets(mapping_traits_genes.traits_to_genes, phegeni_genes, genes_to_proteins, phegeni_proteins, adjacency_list_proteins);

		const auto [proteins_to_proteoforms, proteoforms_to_proteins] = loadMappingProteinsProteoforms(path_file_proteoform_search);
		const auto phegeni_proteoforms = deductProteoformsFromProteins(proteins_to_proteoforms, phegeni_proteins);
		const um<string, bitset<PHEGENI_PROTEOFORMS>> traits_to_proteoforms = convertProteinSets(traits_to_proteins, phegeni_proteins, proteins_to_proteoforms, phegeni_proteoforms, adjacency_list_proteoforms);

		// Calculate overlaps
		const auto overlapping_gene_set_pairs = findOverlappingPairs(mapping_traits_genes.traits_to_genes, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE);
		cout << "Calculating protein sets overlap..." << endl;
		const auto overlapping_protein_set_pairs = findOverlappingPairs(traits_to_proteins, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE, MIN_SET_SIZE, MAX_SET_SIZE);
		cout << "Calculating proteoform sets overlap..." << endl;
		const auto overlapping_proteoform_set_pairs = findOverlappingPairs(traits_to_proteoforms, 1, PHEGENI_PROTEOFORMS, 1, PHEGENI_PROTEOFORMS);
		const auto non_overlapping_proteoform_set_pairs = findOverlappingPairs(traits_to_proteoforms, 0, 0, 1, PHEGENI_PROTEOFORMS);

		cerr << "Total proteoforms sets: " << traits_to_proteoforms.size() << "\n";
		cerr << "Parejas: " << (traits_to_proteoforms.size() * traits_to_proteoforms.size() - traits_to_proteoforms.size()) / 2 << "\n";
		cerr << "Overlapping: " << overlapping_proteoform_set_pairs.size() << "\n";
		cerr << "Non-overlapping: " << non_overlapping_proteoform_set_pairs.size() << "\n";
		cerr << "Total pairs defined: " << overlapping_proteoform_set_pairs.size() + non_overlapping_proteoform_set_pairs.size() << "\n";

		cout << "Finding examples of gene level only overlap...\n";
		const auto examples = findGeneLevelOnlyOverlapPairs(traits_to_proteoforms,
			overlapping_gene_set_pairs,
			overlapping_protein_set_pairs,
			non_overlapping_proteoform_set_pairs);

		cout << "Writing report...\n";
		ofstream report(path_file_report_trait.data());

		writePhenotypeReport(report, examples, traits_to_names, mapping_traits_genes.traits_to_genes, traits_to_proteins, traits_to_proteoforms,
			phegeni_genes, phegeni_proteins, phegeni_proteoforms);
	}

	// Find pairs of modules/pathways that overlap on gene or protein network, but not in proteoform network
	void doAnalysis(string_view path_file_gene_search,
		string_view path_file_protein_search,
		string_view path_file_proteoform_search,
		string_view path_file_PheGenI,
		string_view path_file_mapping_proteins_to_genes,
		string_view path_file_report_pathway,
		string_view path_file_report_trait) {
		cout << "Searching for gene level only overlap examples..." << endl;

		reportPathwayPairs(path_file_gene_search,
			path_file_protein_search,
			path_file_proteoform_search,
			path_file_report_pathway);
		cout << "\n\n************************************\n\n";
		reportPhenotypePairs(path_file_gene_search,
			path_file_protein_search,
			path_file_proteoform_search,
			path_file_PheGenI,
			path_file_mapping_proteins_to_genes,
			path_file_report_trait);
	}

}  // namespace rule_out_gene_centric_overlap
