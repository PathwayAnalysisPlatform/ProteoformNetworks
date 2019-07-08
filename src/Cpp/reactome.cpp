#include "reactome.hpp"

using namespace std;

// Read list of entities from a PathwayMatcher search file
vs loadReactomeEntitiesIndexToEntitites(std::string_view path_search_file) {
	vs result;
	std::ifstream file_search(path_search_file.data());
	std::string field, entity;
	uss temp_set;

	getline(file_search, field);					// Skip header line
	while (getline(file_search, entity, '\t')) {	// Read Entity (Gene | Protein | Proteoform)
		getline(file_search, field);				// Read rest of line
		temp_set.insert(entity);
	}
	return convert_uss_to_vs(temp_set);
}

load_mapping_proteins_proteoforms_result loadMappingProteinsProteoforms(std::string_view path_file_mapping) {
	ummss protein_to_proteforms;
	ummss proteoform_to_proteins;
	std::ifstream file_mapping(path_file_mapping.data());
	std::string protein, proteoform, leftover;

	if (!file_mapping.is_open()) {
		throw std::runtime_error("Error reading file at: " __FUNCTION__);
	}

	getline(file_mapping, protein);  // Discard the header line.
	while (getline(file_mapping, proteoform, '\t')) {
		getline(file_mapping, protein, '\t');
		getline(file_mapping, leftover);
		protein_to_proteforms.emplace(protein, proteoform);
		proteoform_to_proteins.emplace(proteoform, protein);
	}

	std::cerr << "Mapping proteins proteoforms loaded: " << protein_to_proteforms.size() << " = " << proteoform_to_proteins.size() << "\n";
	return { protein_to_proteforms, proteoform_to_proteins };
}

// Create bimap from entities in a PathwayMatcher search file
bimap loadReactomeEntities(std::string_view path_search_file) {
	auto index_to_entities = loadReactomeEntitiesIndexToEntitites(path_search_file.data());
	return createBimap(index_to_entities);
}

reactome_gene_sets loadGeneSets(std::string_view file_path, const bimap& genes, bool pathways) {
	reactome_gene_sets result;
	std::ifstream file_search(file_path.data());
	std::string field, gene, reaction, pathway;

	getline(file_search, field);                // Skip csv header line
	while (getline(file_search, gene, '\t')) {  // Read the members of each gene_set // Read GENE
		getline(file_search, field, '\t');       // Read UNIPROT
		getline(file_search, reaction, '\t');    // Read REACTION_STID
		getline(file_search, field, '\t');       // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');     // Read PATHWAY_STID
		getline(file_search, field);             // Read PATHWAY_DISPLAY_NAME

		result[pathways ? pathway : reaction].set(genes.indexes.find(gene)->second);
	}
	return result;
}

reactome_protein_sets loadProteinSets(std::string_view file_path, const bimap& proteins, bool pathways) {
	reactome_protein_sets result;
	std::ifstream file_search(file_path.data());
	std::string field, entity, reaction, pathway;

	getline(file_search, field);                  // Skip csv header line
	while (getline(file_search, entity, '\t')) {  // Read the members of each gene_set // Read UNIPROT
		getline(file_search, reaction, '\t');      // Read REACTION_STID
		getline(file_search, field, '\t');         // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');       // Read PATHWAY_STID
		getline(file_search, field);               // Read PATHWAY_DISPLAY_NAME

		result[pathways ? pathway : reaction].set(proteins.indexes.find(entity)->second);
	}
	return result;
}

reactome_proteoform_sets loadProteoformSets(std::string_view file_path, const bimap& proteoforms, bool pathways) {
	reactome_proteoform_sets result;
	std::ifstream file_search(file_path.data());
	std::string field, entity, reaction, pathway;

	getline(file_search, field);                  // Skip csv header line
	while (getline(file_search, entity, '\t')) {  // Read the members of each gene_set // Read PROTEOFORM
		getline(file_search, field, '\t');         // Read UNIPROT
		getline(file_search, reaction, '\t');      // Read REACTION_STID
		getline(file_search, field, '\t');         // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');       // Read PATHWAY_STID
		getline(file_search, field);               // Read PATHWAY_DISPLAY_NAME

		result[pathways ? pathway : reaction].set(proteoforms.indexes.find(entity)->second);
	}
	return result;
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

bimap deductProteinsFromGenes(
	std::string_view path_file_mapping_proteins_genes,
	const ummss& genes_to_proteins,
	const bimap& genes) {
	uss temp_protein_set;

	for (const auto& gene_entry : genes.indexes) {
		auto ret = genes_to_proteins.equal_range(gene_entry.first);
		for (auto it = ret.first; it != ret.second; ++it) {
			temp_protein_set.insert(it->second);
		}
	}

	vs index_to_proteins = convert_uss_to_vs(temp_protein_set);
	bimap proteins = createBimap(index_to_proteins);

	cerr << "PHEGEN proteins: " << proteins.entities.size() << " = " << proteins.indexes.size() << "\n";

	return proteins;
}

bimap deductProteoformsFromProteins(const ummss& proteins_to_proteoforms, const bimap& proteins) {
	uss temp_proteoform_set;

	for (const auto& entry : proteins.indexes) {
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

	auto index_to_proteoforms = convert_uss_to_vs(temp_proteoform_set);
	bimap proteoforms = createBimap(index_to_proteoforms);

	cerr << "PHEGEN proteoforms: " << proteoforms.entities.size() << " = " << proteoforms.indexes.size() << "\n";

	return proteoforms;
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
	const auto proteins = createBimap(search_file_path);
	const auto reactions_to_proteins = loadProteinSets(search_file_path, proteins, false);

	for (const auto& reaction_entry : reactions_to_proteins) {
		vs members;
		for (int I = 0; I < reaction_entry.second.size(); I++) {
			if (reaction_entry.second.test(I)) {
				members.push_back(proteins.entities[I]);
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
	const auto proteoforms = createBimap(search_file_path);
	const auto reactions_to_proteoforms = loadProteoformSets(search_file_path, proteoforms, false);

	for (const auto& reaction_entry : reactions_to_proteoforms) {
		vs members;
		for (int I = 0; I < reaction_entry.second.size(); I++) {
			if (reaction_entry.second.test(I)) {
				members.push_back(proteoforms.entities[I]);
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