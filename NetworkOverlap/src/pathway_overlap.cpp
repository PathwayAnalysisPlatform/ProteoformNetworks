//============================================================================
// Description : Get pathway pairs with different overlap for genes, proteins
//               and proteoforms.
// Author      : Luis Francisco Hernández Sánchez
// Copyright   : Licence Apache 2.0
//============================================================================

#include "pathway_overlap.h"

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
const string PATH_FILE_RATIOS = PATH_PHENI + "modified_percentage.csv";

void CreateRatiosFile() {

	ifstream trait_file(PATH_PHENI + "phenotypes.txt");
	ofstream file_modified_percentage(PATH_FILE_RATIOS);

	string proteoform;

	std::regex expression { "[;,]\\d{5}" };

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

					std::smatch modification;
					if (std::regex_search(proteoform, modification, expression)) {
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

set<string> GetVertices(Entity entity, string pathway) {

	set<string> vertices;
	string vertex_name;															// Entity label: protein accession, gene name, etc.
	string line_leftover;
	string path_file_vertices = PATH_DISEASE_MODULES + pathway + "/" + entity_str[entity] + "Vertices.tsv";
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

set<string> GetOverlap(Entity entity, const string& one_trait, const string& other_trait) {

//	cout << "getOverlap: " << entity_str[entity] << "\t" << one_trait << "\t" << other_trait << endl;

	set<string> one_trait_vertices = GetVertices(entity, one_trait);
	set<string> other_trait_vertices = GetVertices(entity, other_trait);
	set<string> overlap;
	set_intersection(one_trait_vertices.begin(),
								one_trait_vertices.end(),
								other_trait_vertices.begin(),
								other_trait_vertices.end(),
								std::inserter(overlap, overlap.begin()));
	return overlap;
}

void CreatePairsFile(float min_modified_percentage) {

	ifstream file_modified_percentage(PATH_FILE_RATIOS);
	if(!file_modified_percentage.good()) {
		CreateRatiosFile();
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
					<< "Difference" << endl;

	for(const auto &one_trait : trait_to_ratio_modified) {
		if (one_trait.second >= min_modified_percentage) {						// For all phenotypes that have percentage higher than threshold
			for(const auto &other_trait : trait_to_ratio_modified) { 			// Compare to all other phenotypes for overlap
				if(one_trait == other_trait) {
					continue;
				}
				cout << "Comparing combination " << ++combination << ": " << one_trait.first << " with " << other_trait.first << endl;

				// If overlap in proteoforms is different than in proteins, add to candidate list
				set<string> overlap_proteins = GetOverlap(protein, one_trait.first, other_trait.first);
				set<string> overlap_proteoforms = GetOverlap(proteoform, one_trait.first, other_trait.first);

				if(overlap_proteins.size() != overlap_proteoforms.size()) {
					unsigned diff = abs((int)overlap_proteins.size() - (int)overlap_proteoforms.size());
					cout << "-- Found candidate: " << one_trait.first << " and " << other_trait.first << " with difference " << diff << endl;
					pairs_diff_overlap.insert(std::make_pair(one_trait.first, other_trait.first));
					file_pairs_diff_overlap << one_trait.first << "\t"
									<< other_trait.first << "\t"
									<< one_trait.second << "\t"
									<< other_trait.second << "\t"
									<< overlap_proteins.size() << "\t"
									<< overlap_proteoforms.size() << "\t"
									<< diff << endl;
				}
			}
		}
	}

	cout << " --------------- pairs_diff_overlap -------------------- " << endl << endl;

	for(const auto &candidate : pairs_diff_overlap) {
		cout << candidate.first << " -- " << candidate.second << endl;
	}
}
