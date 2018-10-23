//============================================================================
// Description : Get phenotype pairs that overlap only with modified proteins,
// 				 or that artefactually overlap at the gene or protein level
// Author      : Luis Francisco Hernández Sánchez
// Copyright   : Licence Apache 2.0
//============================================================================

#include "phenotype_overlap.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;

const string PATH_PHENI = "../resources/PheGenI/";
const string PATH_DISEASE_MODULES = "../diseaseModules/";

const std::regex modification_rgx { "[;,]\\d{5}" };

bool IsModified(const string& proteoform) {
	std::smatch modification;
	return std::regex_search(proteoform, modification, modification_rgx);
}

set<string> GetTraits() {
	ifstream trait_file(PATH_PHENI + "traits.txt");
	set<string> traits;

	if(trait_file) {
		string trait = "";
		getline(trait_file, trait); 						// Discard the first line of header
		while (getline(trait_file, trait)) {
			traits.insert(trait);
		}
	}

	return traits;
}

void CreatePhenotypeRatiosFile() {

	ifstream trait_file(PATH_PHENI + "Phenotype.txt");
	ofstream file_modified_percentage(PATH_PHENI + "modified_percentage.csv");

	string proteoform;

	file_modified_percentage << "Trait" << "\t"
							 << "Num_proteoforms" << "\t"
							 << "Num_modified" << "\t"
							 << "Ratio" << endl;

	if (trait_file.good()) {
		string trait = "";
		getline(trait_file, trait); 						// Discard the first line of header
		while (getline(trait_file, trait)) {

			int num_modified = 0;
			int num_proteoforms = 0;
			double ratio = 0.0;
			ifstream file_vertices_proteoforms(PATH_DISEASE_MODULES + trait + "/proteoformVertices.tsv");

			if (file_vertices_proteoforms.good()) {
				getline(file_vertices_proteoforms, proteoform);	// Skip the first line
				while (getline(file_vertices_proteoforms, proteoform, '\t')) {
					num_proteoforms++;

					string rest_of_line;
					getline(file_vertices_proteoforms, rest_of_line);

					if (IsModified(proteoform)) {
						num_modified++;
					}
				}
			} else {
				cout << "Failed to open the proteoforms file." << endl;
			}

			ratio = num_proteoforms ? (double) num_modified / (double) num_proteoforms : 0.0;

			cout << trait << "\t"
				 << num_proteoforms << "\t"
				 << num_modified << "\t"
				 << ratio << endl;
			file_modified_percentage << trait << "\t"
									 << num_proteoforms << "\t"
									 << num_modified << "\t"
									 << ratio << endl;
		}
	}
}

set<string> GetPhenotypeVertices(Entity entity, const string& trait) {

	set<string> vertices;
	string vertex_name;															// Entity label: protein accession, gene name, etc.
	string line_leftover;
	string path_file_vertices = PATH_DISEASE_MODULES + trait + "/" + entity_str[entity] + "Vertices.tsv";
	ifstream file_vertices(path_file_vertices);

	if(file_vertices) {
		getline(file_vertices, line_leftover); 									// Discard the first line
		while(file_vertices >> vertex_name) {
			getline(file_vertices, line_leftover);
//			cout << " ++ " << vertex_name << endl;
			vertices.insert(vertex_name);
		}
	}
//	cout << "++++++++++" << endl;
	return vertices;
}

set<string> GetPhenotypeVerticesFromEdgesFile(Entity entity, const string& trait) {
	set<string> vertices;

	string one_vertex, other_vertex;
	string line_leftover;
	string path_file_edges = PATH_DISEASE_MODULES + trait + "/" + entity_str[entity] + "InternalEdges.tsv";
	ifstream file_edges(path_file_edges);

	if(file_edges) {
		getline(file_edges, line_leftover);										//Discard first line with the header
		while(file_edges >> one_vertex >> other_vertex) {
			getline(file_edges, line_leftover);
			vertices.insert(one_vertex);
			vertices.insert(other_vertex);
		}
	}

	return vertices;
}

multimap<string, string> GetPhenotypeEdges(Entity entity, const string& trait) {
	multimap<string, string> edges;

	string one_vertex, other_vertex;
	string line_leftover;
	string path_file_edges = PATH_DISEASE_MODULES + trait + "/" + entity_str[entity] + "InternalEdges.tsv";
	ifstream file_edges(path_file_edges);

	if(file_edges) {
		getline(file_edges, line_leftover);										//Discard first line with the header
		while(file_edges >> one_vertex >> other_vertex) {
			getline(file_edges, line_leftover);
			edges.insert(make_pair(one_vertex, other_vertex));
			edges.insert(make_pair(other_vertex, one_vertex));
		}
	}

	return edges;
}

set<string> GetVertexOverlap(Entity entity, const string& one_trait, const string& other_trait) {

	set<string> one_trait_vertices, other_trait_vertices, overlap;

	one_trait_vertices = GetPhenotypeVerticesFromEdgesFile(entity, one_trait);	// Reads from edges file, to keep only connected vertices
	other_trait_vertices = GetPhenotypeVerticesFromEdgesFile(entity, other_trait);

	set_intersection(one_trait_vertices.begin(),
					 one_trait_vertices.end(),
					 other_trait_vertices.begin(),
					 other_trait_vertices.end(),
					 std::inserter(overlap, overlap.begin()));

	return overlap;
}

set<string> GetVertexOverlap(Entity entity, const string& one_trait, const string& other_trait, bool verbose) {

	set<string> overlap = GetVertexOverlap(entity, one_trait, other_trait);

	if(verbose) {
		cout << " ----- " << entity_str[entity] << " ----- " << endl;
		for(const auto &entity : overlap) {
			cout << entity << endl;
		}
		cout << " *********************** " << endl;
	}

	return overlap;
}

set<string> GetVertexOverlap(Entity entity, set<string> one_trait_vertices, set<string> other_trait_vertices) {

	set<string> overlap;

	set_intersection(one_trait_vertices.begin(),
					 one_trait_vertices.end(),
					 other_trait_vertices.begin(),
					 other_trait_vertices.end(),
					 std::inserter(overlap, overlap.begin()));

	return overlap;
}

void CreatePhenotypePairsFile(double min_modified_percentage, Entity entity_type) {

	// Find the phenotype pairs of two types:
	//		* Type 1: Overlapping nodes are only modified proteoforms
	// 		* Type 2: They overlap in the gene or protein network, but not in the proteoform network

	ifstream file_modified_percentage(PATH_PHENI + "modified_percentage.csv");
	if(!file_modified_percentage.good()) {
		CreatePhenotypeRatiosFile();
	}

	char trait_arr[80];
	int num_proteoforms = -1;
	int num_modified = -1;
	int combination = 0;
	float ratio = -1.0;
	string trait;

	multimap<string, string> pairs_diff_overlap;										// Disease pairs with different overlap size for proteins and proteoform networks
	map<string, double> trait_to_ratio_modified; 											// Percentage of modified proteins

	freopen((PATH_PHENI + "modified_percentage.csv").c_str(), "r", stdin);
	ofstream file_pairs_diff_overlap(PATH_PHENI + "pairs_diff_overlap.csv");

	// Read all traits
	getline(cin, trait); // Discard file header
	while(scanf("%s%i%d%f", trait_arr, &num_proteoforms, &num_modified, &ratio) > 1) {
//		cout << trait_arr << " ** " << num_proteoforms << " ** " << num_modified << " ** " << ratio << endl;
		trait.assign(trait_arr);
		trait_to_ratio_modified[trait] = ratio;
	}

	fclose(stdin);

	file_pairs_diff_overlap << "Trait1\t"
					<< "Trait2\t"
					<< "RatioModified1\t"
					<< "RatioModified2\t"
					<< "OverlapProteins\t"
					<< "OverlapProteoforms\t"
					<< "Type" << endl;

	for(const auto &one_trait : trait_to_ratio_modified) {
		if (one_trait.second >= min_modified_percentage) {						// For all Phenotype that have percentage higher than threshold

			set<string> one_trait_entity_type_vertices = GetPhenotypeVertices(entity_type, one_trait.first);	// Reads from edges file, to keep only connected vertices
			set<string> one_trait_proteoform_vertices = GetPhenotypeVerticesFromEdgesFile(proteoform, one_trait.first);

			for(const auto &other_trait : trait_to_ratio_modified) { 			// Compare to all other Phenotype for overlap
				if(one_trait == other_trait) {									// Discard pairs a -- a
					continue;
				}

				if(other_trait.second >= min_modified_percentage) {				// Discard inverted pairs to avoid, having a -- b and b -- a
					if(one_trait.first > other_trait.first) {
						continue;
					}
				}
				cout << "Comparing " << ++combination << ": " << one_trait.first << " -- " << other_trait.first << endl;

				set<string> other_trait_entity_type_vertices = GetPhenotypeVertices(entity_type, other_trait.first);
				set<string> other_trait_proteoform_vertices = GetPhenotypeVerticesFromEdgesFile(proteoform, other_trait.first);

				// If overlap in proteoforms is different than in proteins, add to candidate list
				set<string> overlap_entity_type = GetVertexOverlap(entity_type, one_trait_entity_type_vertices, other_trait_entity_type_vertices);
				set<string> overlap_proteoforms = GetVertexOverlap(proteoform, one_trait_proteoform_vertices, other_trait_proteoform_vertices);

				if(overlap_entity_type.size() != overlap_proteoforms.size()) {

					// Check for type 1: If all overlapping proteoforms are modified
					bool overlaps_only_with_modified_proteoforms = overlap_proteoforms.size();
					for(const auto &proteoform : overlap_proteoforms) {
						if(!IsModified(proteoform)) {
							overlaps_only_with_modified_proteoforms = false;
						}
					}

					if(overlaps_only_with_modified_proteoforms) {
						cout << "-- Found candidate type 1: " << one_trait.first << " and " << other_trait.first << " opo" << endl;
						pairs_diff_overlap.insert(std::make_pair(one_trait.first, other_trait.first));
						file_pairs_diff_overlap << one_trait.first << "\t"
												<< other_trait.first << "\t"
												<< one_trait.second << "\t"
												<< other_trait.second << "\t"
												<< overlap_entity_type.size() << "\t"
												<< overlap_proteoforms.size() << "\t"
												<< "opo" << endl;
					}

					// Check for type 2: If diseases do not overlap in the proteoform network but overlap in the other network
					if(overlap_entity_type.size() && !overlap_proteoforms.size()) {
						cout << "-- Found candidate: " << one_trait.first << " and " << other_trait.first << " sep" << endl;
						pairs_diff_overlap.insert(std::make_pair(one_trait.first, other_trait.first));
						file_pairs_diff_overlap << one_trait.first << "\t"
										<< other_trait.first << "\t"
										<< one_trait.second << "\t"
										<< other_trait.second << "\t"
										<< overlap_entity_type.size() << "\t"
										<< overlap_proteoforms.size() << "\t"
										<< "sep" << endl;
					}
				}
			}
		}
	}

	cout << " --------------- pairs_diff_overlap -------------------- " << endl << endl;

	for(const auto &candidate : pairs_diff_overlap) {
		cout << candidate.first << " -- " << candidate.second << endl;
	}
}

void CheckNumVerticesForTraits() {
	// Check number of connected vertices and total number vertices in the graph

	std::set<std::string> traits = GetTraits();

	for(const auto &trait : traits) {

		std::set<std::string> vertices_from_vertices_file = GetPhenotypeVertices(gene, trait);
		std::set<std::string> vertices_from_edges_file = GetPhenotypeVerticesFromEdgesFile(gene, trait);

		if(vertices_from_vertices_file.size() != vertices_from_edges_file.size()) {
			std::cout << "Difference in " << trait << "\t\t\t";
			std::cout << "from vertices: " << vertices_from_vertices_file.size();
			std::cout << "\tfrom edges: " << vertices_from_edges_file.size() << std::endl;
		}
	}
}

multimap<string, string> GetPairsOverlappingAtEntities(set<string> entities, Entity entity_type, multimap<string, string> candidate_pairs) {

	multimap<string, string> pairs;

	for(const auto &pair : candidate_pairs) {
		set<string> overlap = GetVertexOverlap(entity_type, pair.first, pair.second);

			bool contains_all = true;
			for(const auto &entity : entities) {
				if(overlap.find(entity) == overlap.end()) {
					contains_all = false;
					break;
				}
			}

			if(contains_all) {
				pairs.insert(std::make_pair(pair.first, pair.second));
			}
	}

	return pairs;
}

set<string> GetProteoformsWithAccession(const string& protein_accession, set<string> candidate_proteoforms) {
	set<string> proteoforms;
	for(const auto &proteoform : candidate_proteoforms) {
		if(proteoform.find(protein_accession) != string::npos) {
			proteoforms.insert(proteoform);
		}
	}
	return proteoforms;
}

void ShowOverlapsAndSplittingProteoforms(std::multimap<std::string, std::string> pairs) {
	// For each pair of entity sets, show why they are not overlapping in the proteoform network

	for(const auto &pair : pairs) {
		string one_trait = pair.first;
		string other_trait = pair.second;

		cout << " ---- " << one_trait << " ---- " << other_trait << " ---- " << endl;

		set<string> overlap_genes = GetVertexOverlap(gene, one_trait, other_trait, true);
		set<string> overlap_proteins = GetVertexOverlap(protein, one_trait, other_trait, true);
		set<string> overlap_proteoforms = GetVertexOverlap(proteoform, one_trait, other_trait, true);

		// Show the members of the proteoform network with that accession
		set<string> proteoforms_splitting;

		cout << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
		cout << one_trait << " has: " << endl;
		for(const auto &protein : overlap_proteins) {
			set<string> proteoforms = GetProteoformsWithAccession(protein, GetPhenotypeVerticesFromEdgesFile(proteoform, one_trait));
			for(const auto &proteoform : proteoforms) {
				cout << proteoform << endl;
			}
		}
		cout << " <<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
		cout << other_trait << " has: " << endl;
		for(const auto &protein : overlap_proteins) {
			set<string> proteoforms = GetProteoformsWithAccession(protein, GetPhenotypeVerticesFromEdgesFile(proteoform, other_trait));
			for(const auto &proteoform : proteoforms) {
				cout << proteoform << endl;
			}
		}
		cout << " <<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>" << endl << endl;
	}
}
