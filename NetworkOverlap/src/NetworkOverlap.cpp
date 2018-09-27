//============================================================================
// Name        : NetworkOverlap.cpp
// Author      : Luis Francisco Hernández Sánchez
// Version     : 0.0.2
// Copyright   : Licence Apache 2.0
// Description : Get trait pairs with different network overlap for proteins
//               and proteoforms.
//============================================================================

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;

enum Entity { gene, protein, proteoform };
vector<string> entity_str = {"gene", "protein", "proteoform"};

const float MIN_PERCENTAGE_MODIFIED {0.9};
const string PATH_PHENI = "../resources/PheGenI/";
const string PATH_DISEASE_MODULES = "../diseaseModules/";

void create_ratios_file() {

	ifstream trait_file((PATH_PHENI + "traits.txt").c_str());;
	ofstream file_modified_percentage((PATH_PHENI + "modified_percentage.csv").c_str());

	string gene = "";
	string proteoform = "";

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
			ifstream file_vertices_proteoforms((PATH_DISEASE_MODULES + trait + "/proteoformVertices.tsv").c_str());

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

set<string> getVertices(Entity entity, string trait) {

	set<string> vertices;
	string vertex_name;															// Entity label: protein accession, gene name, etc.
	string line_leftover;
	string path_file_vertices = PATH_DISEASE_MODULES + trait + "/" + entity_str[entity] + "Vertices.tsv";
	ifstream file_vertices(path_file_vertices.c_str());

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

set<string> getOverlap(Entity entity, const string& one_trait, const string& other_trait) {

//	cout << "getOverlap: " << entity_str[entity] << "\t" << one_trait << "\t" << other_trait << endl;

	set<string> one_trait_vertices = getVertices(entity, one_trait);
	set<string> other_trait_vertices = getVertices(entity, other_trait);
	set<string> overlap;
	set_intersection(one_trait_vertices.begin(),
								one_trait_vertices.end(),
								other_trait_vertices.begin(),
								other_trait_vertices.end(),
								std::inserter(overlap, overlap.begin()));
	return overlap;
}

int main() {

//	create_ratios_file();

	char trait_arr[80];
	int num_proteoforms = -1;
	int num_modified = -1;
	int combination = 0;
	float ratio = -1.0;
	string trait{};

	multimap<string, string> candidates;										// Disease pairs with different overlap size for proteins and proteoform networks
	map<string, double> trait_to_ratio_modified; 								// Percentage of modified proteins

	freopen((PATH_PHENI + "modified_percentage.csv").c_str(), "r", stdin);
	ofstream file_candidates(PATH_PHENI + "candidate_disease_pairs.csv");

	// Read all traits
	getline(cin, trait); // Discard file header
	while(scanf("%s%i%d%f", trait_arr, &num_proteoforms, &num_modified, &ratio) > 1) {
//		cout << trait_arr << " ** " << num_proteoforms << " ** " << num_modified << " ** " << ratio << endl;
		trait.assign(trait_arr);
		trait_to_ratio_modified[trait] = ratio;
	}

	fclose(stdin);

	file_candidates << "Trait1\t"
					<< "Trait2\t"
					<< "RatioModified1\t"
					<< "RatioModified2\t"
					<< "OverlapProteins\t"
					<< "OverlapProteoforms\t"
					<< "Difference" << endl;

	for(const auto &one_trait : trait_to_ratio_modified) {
		if (one_trait.second >= MIN_PERCENTAGE_MODIFIED) {						// For all phenotypes that have percentage higher than threshold
			for(const auto &other_trait : trait_to_ratio_modified) { 			// Compare to all other phenotypes for overlap
				if(one_trait == other_trait) {
					continue;
				}
				cout << "Comparing combination " << ++combination << ": " << one_trait.first << " with " << other_trait.first << endl;

				// If overlap in proteoforms is different than in proteins, add to candidate list
				set<string> overlap_proteins = getOverlap(protein, one_trait.first, other_trait.first);
				set<string> overlap_proteoforms = getOverlap(proteoform, one_trait.first, other_trait.first);

				if(overlap_proteins.size() != overlap_proteoforms.size()) {
					unsigned diff = abs((int)overlap_proteins.size() - (int)overlap_proteoforms.size());
					cout << "-- Found candidate: " << one_trait.first << " and " << other_trait.first << " with difference " << diff << endl;
					candidates.insert(std::make_pair(one_trait.first, other_trait.first));
					file_candidates << one_trait.first << "\t"
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

	cout << " --------------- Candidates -------------------- " << endl << endl;

	for(const auto &candidate : candidates) {
		cout << candidate.first << " -- " << candidate.second << endl;
	}
	return 0;
}
