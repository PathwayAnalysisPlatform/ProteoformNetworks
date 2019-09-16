#include "reactome.hpp"

using namespace std;

umsb loadGeneSets(string_view file_path, const bimap_str_int& phegeni_genes, bool pathways) {
    umsb result;
    // TODO: Initialize result structure to the bitsets
	ifstream file_search(file_path.data());
	string field, gene, reaction, pathway;

	getline(file_search, field);                // Skip csv header line
	while (getline(file_search, gene, '\t')) {  // Read the members of each gene_set // Read GENE
		getline(file_search, field, '\t');       // Read UNIPROT
		getline(file_search, reaction, '\t');    // Read REACTION_STID
		getline(file_search, field, '\t');       // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');     // Read PATHWAY_STID
		getline(file_search, field);             // Read PATHWAY_DISPLAY_NAME

		result[pathways ? pathway : reaction][phegeni_genes.str_to_int.find(gene)->second].set();
	}
	return result;
}



umss loadPathwayNames(string_view path_search_file) {
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

bimap_str_int deductProteinsFromGenes(
	const ummss& genes_to_proteins,
	const bimap_str_int& phegeni_genes) {
	uss temp_protein_set;

	for (const auto& gene_entry : phegeni_genes.str_to_int) {
		auto ret = genes_to_proteins.equal_range(gene_entry.first);
		for (auto it = ret.first; it != ret.second; ++it) {
			temp_protein_set.insert(it->second);
		}
	}

	vs index_to_proteins = convert_uss_to_vs(temp_protein_set);
	bimap_str_int proteins = createBimap(index_to_proteins);

	std::cerr << "PHEGEN proteins: " << proteins.int_to_str.size() << " = " << proteins.str_to_int.size() << "\n";

	return proteins;
}

bimap_str_int deductProteoformsFromProteins(const ummss& proteins_to_proteoforms, const bimap_str_int& proteins) {
	uss temp_proteoform_set;

	for (const auto& entry : proteins.str_to_int) {
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
	bimap_str_int proteoforms = createBimap(index_to_proteoforms);

	cerr << "PHEGEN proteoforms: " << proteoforms.int_to_str.size() << " = " << proteoforms.str_to_int.size() << "\n";

	return proteoforms;
}

// Read entity interactions in Reactome from a PathwayMatcher edges file
ummss loadInteractions(std::string_view path_file_edges) {

    // TODO
	ummss adjacenty_list;
//	const auto proteoforms = createBimap(search_file_path);
//	const auto reactions_to_proteoforms = loadProteoformSets(search_file_path, proteoforms, false);
//
//	for (const auto& reaction_entry : reactions_to_proteoforms) {
//		vs members;
//		for (int I = 0; I < reaction_entry.second.size(); I++) {
//			if (reaction_entry.second.test(I)) {
//				members.push_back(proteoforms.int_to_str[I]);
//			}
//		}
//
//		for (const auto& one_member : members) {
//			for (const auto& other_member : members) {
//				adjacenty_list.emplace(one_member, other_member);
//			}
//		}
//	}

	return adjacenty_list;
}