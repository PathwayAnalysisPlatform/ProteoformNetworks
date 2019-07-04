#include "overlap.hpp"

using namespace std;

void doOverlapAnalysis(std::string_view path_file_PheGenI, std::string_view path_file_genes) {

	cout << "Load Reactome gene set\n";
	const bimap reactome_genes = createBimap(path_file_genes, true);

	cout << "Loading PheGen gene sets\n";
	const auto [phegeni_genes, phegeni_traits] = loadPheGenIEntities(path_file_PheGenI, GENOME_WIDE_SIGNIFICANCE, reactome_genes.indexes);

	std::cout << "Number of genes in Reactome: " << reactome_genes.entities.size() << "\n";
	std::cout << "Number of genes in PheGenI: " << phegeni_genes.entities.size() << "\n";
	std::cout << "Number of traits in PheGenI: " << phegeni_traits.entities.size() << "\n";
}

unordered_map<string, bitset<NUM_GENES>> loadGeneSets(std::string_view file_path, const umsi& entities_to_index, bool pathways) {
	unordered_map<string, bitset<NUM_GENES>> result;
	ifstream file_search(file_path.data());
	string field, gene, reaction, pathway;

	getline(file_search, field);                // Skip csv header line
	while (getline(file_search, gene, '\t')) {  // Read the members of each gene_set // Read GENE
		getline(file_search, field, '\t');       // Read UNIPROT
		getline(file_search, reaction, '\t');    // Read REACTION_STID
		getline(file_search, field, '\t');       // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');     // Read PATHWAY_STID
		getline(file_search, field);             // Read PATHWAY_DISPLAY_NAME

		result[pathways ? pathway : reaction].set(entities_to_index.find(gene)->second);
	}
	return result;
}

umss createTraitNames(const unordered_map<string, bitset<NUM_PHEGEN_GENES>>& traits_to_genes) {
	umss sets_to_names;
	for (const auto& trait_entry : traits_to_genes) {
		sets_to_names.emplace(trait_entry.first, trait_entry.first);
	}
	return sets_to_names;
}

unordered_map<string, bitset<NUM_PROTEINS>> loadProteinSets(std::string_view file_path, const umsi& entities_to_index, bool pathways) {
	unordered_map<string, bitset<NUM_PROTEINS>> result;
	ifstream file_search(file_path.data());
	string field, entity, reaction, pathway;

	getline(file_search, field);                  // Skip csv header line
	while (getline(file_search, entity, '\t')) {  // Read the members of each gene_set // Read UNIPROT
		getline(file_search, reaction, '\t');      // Read REACTION_STID
		getline(file_search, field, '\t');         // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');       // Read PATHWAY_STID
		getline(file_search, field);               // Read PATHWAY_DISPLAY_NAME

		result[pathways ? pathway : reaction].set(entities_to_index.find(entity)->second);
	}
	return result;
}

unordered_map<string, bitset<NUM_PROTEOFORMS>> loadProteoformSets(std::string_view file_path, const umsi& entities_to_index, bool pathways) {
	unordered_map<string, bitset<NUM_PROTEOFORMS>> result;
	ifstream file_search(file_path.data());
	string field, entity, reaction, pathway;

	getline(file_search, field);                  // Skip csv header line
	while (getline(file_search, entity, '\t')) {  // Read the members of each gene_set // Read PROTEOFORM
		getline(file_search, field, '\t');         // Read UNIPROT
		getline(file_search, reaction, '\t');      // Read REACTION_STID
		getline(file_search, field, '\t');         // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');       // Read PATHWAY_STID
		getline(file_search, field);               // Read PATHWAY_DISPLAY_NAME

		result[pathways ? pathway : reaction].set(entities_to_index.find(entity)->second);
	}
	return result;
}

template <size_t total_num_entities>
uss getEntityStrings(const bitset<total_num_entities>& entity_set, const vs& index_to_entities) {
	uss result;
	for (int I = 0; I < total_num_entities; I++) {
		if (entity_set.test(I)) {
			result.insert(index_to_entities[I]);
		}
	}
	return result;
}

uss getGeneStrings(const bitset<NUM_GENES>& gene_set, const vs& index_to_genes) {
	return getEntityStrings(gene_set, index_to_genes);
}

uss getProteinStrings(const bitset<NUM_PROTEINS>& protein_set, const vs& index_to_proteins) {
	return getEntityStrings(protein_set, index_to_proteins);
}

uss getProteoformStrings(const bitset<NUM_PROTEOFORMS>& proteoform_set, const vs& index_to_proteoforms) {
	return getEntityStrings(proteoform_set, index_to_proteoforms);
}

uss getGeneStrings(const bitset<NUM_PHEGEN_GENES>& gene_set, const vs& index_to_genes) {
	return getEntityStrings(gene_set, index_to_genes);
}

uss getProteinStrings(const bitset<NUM_PHEGEN_PROTEINS>& protein_set, const vs& index_to_proteins) {
	return getEntityStrings(protein_set, index_to_proteins);
}

uss getProteoformStrings(const bitset<NUM_PHEGEN_PROTEOFORMS>& proteoform_set, const vs& index_to_proteoforms) {
	return getEntityStrings(proteoform_set, index_to_proteoforms);
}

load_trait_gene_sets_result loadTraitGeneSets(const std::string& path_file_phegen,
	const double& max_p_value,
	const bimap& genes,
	const bimap& traits,
	const umsi& reactome_genes_to_index) {
	ifstream file_phegen(path_file_phegen);
	string line, field, trait, gene, gene2;
	string p_value_str;
	long double p_value;
	um<std::string, std::bitset<NUM_PHEGEN_GENES>> sets_to_genes;
	um<std::string, std::bitset<NUM_PHEGEN_TRAITS>> genes_to_sets;

	if (!file_phegen.is_open()) {
		throw runtime_error("Cannot open " + path_file_phegen);
	}

	getline(file_phegen, line);                  // Read header line
	while (getline(file_phegen, field, '\t')) {  // Read #
		getline(file_phegen, trait, '\t');        // Read Trait
		getline(file_phegen, field, '\t');        // Read SNP rs
		getline(file_phegen, field, '\t');        // Read Context
		getline(file_phegen, gene, '\t');         //	Gene
		getline(file_phegen, field, '\t');        //	Gene ID
		getline(file_phegen, gene2, '\t');        //	Gene 2
		getline(file_phegen, field, '\t');        //	Gene ID 2
		getline(file_phegen, field, '\t');        // Read Chromosome
		getline(file_phegen, field, '\t');        // Read Location
		getline(file_phegen, p_value_str, '\t');  // Read P-Value
		getline(file_phegen, line);               // Skip header line leftoever: Source,	PubMed,	Analysis ID,	Study ID,	Study Name

		try {
			p_value = stold(p_value_str);
			// cerr << "Converted correctly: " << p_value_str << "\n";
			if (p_value <= GENOME_WIDE_SIGNIFICANCE) {
				if (reactome_genes_to_index.find(gene) != reactome_genes_to_index.end()) {
					sets_to_genes[trait].set(genes.indexes.at(gene));
					genes_to_sets[gene].set(traits.indexes.at(trait));
				}
				if (reactome_genes_to_index.find(gene2) != reactome_genes_to_index.end()) {
					sets_to_genes[trait].set(genes.indexes.at(gene2));
					genes_to_sets[gene2].set(traits.indexes.at(trait));
				}
			}
		}
		catch (const std::exception& ex) {
			cerr << "Error converting: **" << p_value_str << "**\n";
		}
	}

	return { genes_to_sets, sets_to_genes };
}

// Loads mapping from one column of strings to a second column of strings. If there are multiple values on the second column for one value on the left, then
//they are separated by spaces. Columns are separated by tabs.
load_mapping_result loadMapping(const std::string& path_file_mapping) {
	ummss ones_to_others;
	ummss others_to_ones;
	ifstream file_mapping(path_file_mapping);
	string to_str, from_str, leftover;

	if (!file_mapping.is_open()) {
		throw runtime_error(path_file_mapping);
	}

	getline(file_mapping, from_str);  // Discard the header line.
	while (getline(file_mapping, from_str, '\t')) {
		char c;

		while (file_mapping.get(c)) {
			if (c == ' ' || c == '\t' || c == '\n') {
				ones_to_others.emplace(from_str, to_str);
				others_to_ones.emplace(to_str, from_str);
				to_str = "";
				if (c == '\t') {
					getline(file_mapping, leftover);
				}
				if (c != ' ') {
					break;
				}
			}
			else {
				to_str = to_str.append(1, c);
			}
		}
	}

	cerr << "Mapping loaded: " << ones_to_others.size() << " = " << others_to_ones.size() << "\n";
	return { ones_to_others, others_to_ones };
}

template <size_t total_num_original_entities, size_t total_num_result_entities>
unordered_map<string, bitset<total_num_result_entities>> convertSets(const unordered_map<string, bitset<total_num_original_entities>>& traits_to_original_entities,
	const vs& index_to_original_entities,
	const ummss& mapping,
	const umsi& result_entities_to_index,
	const ummss& adjacency_list_result_entities) {
	unordered_map<string, bitset<total_num_result_entities>> traits_to_result_entities;
	for (const auto& trait_entry : traits_to_original_entities) {  // For each trait entry
		uss candidates;
		for (int I = 0; I < total_num_original_entities; I++) {  // For each member in this trait
			if (trait_entry.second.test(I)) {
				auto range = mapping.equal_range(index_to_original_entities[I]);  // For each mapping entity
				for (auto it = range.first; it != range.second; it++) {
					candidates.insert(it->second);
				}
			}
		}

		// Keep only those connected to any of the other gene set members in the reference network
		for (const auto& candidate : candidates) {
			auto range = adjacency_list_result_entities.equal_range(candidate);
			for (auto it = range.first; it != range.second; ++it) {
				if (candidates.find(it->second) != candidates.end()) {
					traits_to_result_entities[trait_entry.first].set(result_entities_to_index.at(candidate));  // Set in stone that the candidate is in the new set
					break;
				}
			}
		}
	}

	return traits_to_result_entities;
}

unordered_map<string, bitset<NUM_PHEGEN_PROTEINS>> convertGeneSets(const unordered_map<string, bitset<NUM_PHEGEN_GENES>>& traits_to_genes,
	const bimap& genes,
	const ummss& mapping_genes_to_proteins,
	const umsi& proteins_to_index,
	const ummss& adjacency_list_proteins) {
	cout << "Converting gene to protein sets.\n";
	return convertSets<NUM_PHEGEN_GENES, NUM_PHEGEN_PROTEINS>(traits_to_genes, genes.entities, mapping_genes_to_proteins, proteins_to_index, adjacency_list_proteins);
}

unordered_map<string, bitset<NUM_PHEGEN_PROTEOFORMS>> convertProteinSets(const unordered_map<string, bitset<NUM_PHEGEN_PROTEINS>>& traits_to_proteins,
	const vs& index_to_proteins,
	const ummss& mapping_proteins_to_proteoforms,
	const umsi& proteoforms_to_index,
	const ummss& adjacency_list_proteoforms) {
	cout << "Converting protein to proteoform sets.\n";
	return convertSets<NUM_PHEGEN_PROTEINS, NUM_PHEGEN_PROTEOFORMS>(traits_to_proteins, index_to_proteins, mapping_proteins_to_proteoforms, proteoforms_to_index, adjacency_list_proteoforms);
}

reactome_networks loadReactomeNetworks(std::string_view path_file_gene_search,
	std::string_view path_file_protein_search,
	std::string_view path_file_proteoform_search) {
	cout << "Loading Reactome networks.\n";
	// const ummss adjacency_list_genes = loadGenesAdjacencyList(path_file_gene_search);
	return { loadProteinsAdjacencyList(path_file_protein_search), loadProteoformsAdjacencyList(path_file_proteoform_search) };
}

ummss loadProteinsAdjacencyList(std::string_view search_file_path) {
	ummss adjacenty_list;
	const auto [index_to_entities, entities_to_index] = createBimap(search_file_path);
	const auto reactions_to_entities = loadProteinSets(search_file_path, entities_to_index, false);

	for (const auto& reaction_entry : reactions_to_entities) {
		vs members;
		for (int I = 0; I < reaction_entry.second.size(); I++) {
			if (reaction_entry.second.test(I)) {
				members.push_back(index_to_entities[I]);
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

ummss loadProteoformsAdjacencyList(std::string_view search_file_path) {
	ummss adjacenty_list;
	const auto [index_to_entities, entities_to_index] = createBimap(search_file_path);
	const auto reactions_to_entities = loadProteoformSets(search_file_path, entities_to_index, false);

	for (const auto& reaction_entry : reactions_to_entities) {
		vs members;
		for (int I = 0; I < reaction_entry.second.size(); I++) {
			if (reaction_entry.second.test(I)) {
				members.push_back(index_to_entities[I]);
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

bimap deductProteinsFromGenes(std::string_view path_file_mapping_proteins_genes,
	const umsi& genes_to_index,
	const ummss& genes_to_proteins) {
	uss temp_protein_set;

	for (const auto& gene_entry : genes_to_index) {
		auto ret = genes_to_proteins.equal_range(gene_entry.first);
		for (auto it = ret.first; it != ret.second; ++it) {
			temp_protein_set.insert(it->second);
		}
	}

	vs index_to_proteins = convert_uss_to_vs(temp_protein_set);
	umsi proteins_to_index = createEntitiesToIndex(index_to_proteins);

	cerr << "PHEGEN proteins: " << index_to_proteins.size() << " = " << proteins_to_index.size() << "\n";

	return { index_to_proteins, proteins_to_index };
}



umss loadPathwayNames(std::string_view path_search_file) {
	umss result;
	ifstream file_search(path_search_file.data());
	string field, pathway, pathway_name;
	int field_cont = 0;
	int pathway_stid_column = 0;

	while (getline(file_search, field, '\t')) {
		if (field == "PATHWAY_STID") {
			pathway_stid_column = field_cont;
			break;
		}
		field_cont++;
	}
	getline(file_search, field);  // Skip header line leftoever

	while (getline(file_search, field, '\t')) {                        // While there are records in the file
		for (int column = 1; column < pathway_stid_column; column++) {  // Skip fields before the pathway_stid
			getline(file_search, field, '\t');
		}
		getline(file_search, pathway, '\t');  // Read PATHWAY_STID
		getline(file_search, pathway_name);   // Read PATHWAY_DISPLAY_NAME
		result[pathway] = pathway_name;
	}
	return result;
}

bimap deductProteoformsFromProteins(const ummss& proteins_to_proteoforms, const umsi& proteins_to_index) {
	uss temp_proteoform_set;

	for (const auto& entry : proteins_to_index) {
		auto ret = proteins_to_proteoforms.equal_range(entry.first);
		if (ret.first != ret.second) {
			for (auto it = ret.first; it != ret.second; ++it) {
				temp_proteoform_set.insert(it->second);
			}
		}
		else {
			temp_proteoform_set.insert(entry.first + ";");
		}
	}

	vs index_to_proteoforms = convert_uss_to_vs(temp_proteoform_set);
	umsi proteoforms_to_index = createEntitiesToIndex(index_to_proteoforms);

	cerr << "PHEGEN proteoforms: " << index_to_proteoforms.size() << " = " << proteoforms_to_index.size() << "\n";

	return { index_to_proteoforms, proteoforms_to_index };
}