//============================================================================
// Name        : NetworkOverlap.cpp
// Author      : Luis Francisco Hernández Sánchez
// Version     : 0.0.1
// Copyright   : Licence Apache 2.0
// Description : Calculates the percentage of modified proteins for each trait
//============================================================================

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <regex>
#include <cstdio>
#include <set>
#include <vector>

using namespace std;

enum Entity { gene, protein, proteoform };
vector<string> entity_str = {"gene", "protein", "proteoform"};

const float MIN_PERCENTAGE_MODIFIED {0.25};

string path_phenI = "../resources/PheGenI/";
string path_disease_modules = "../diseaseModules/";
map<string, int> trait_to_num_modified; // Number of modifications for each trait
map<string, double> trait_to_ratio_modified; // Percentage of modified proteins

void create_ratios_file() {
	ifstream trait_file;
	ifstream vertices_file_genes;
	ifstream vertices_file_proteins;

	string gene = "";
	string proteoform = "";

	ofstream file_modified_percentage(
			(path_phenI + "modified_percentage.csv").c_str());

	trait_file.open((path_phenI + "traits.txt").c_str(), std::ifstream::in);

	if (trait_file.good()) {
		string trait = "";
		while (getline(trait_file, trait)) {

			int num_modified = 0;
			int num_proteoforms = 0;
			double ratio = 0.0;
			ifstream file_vertices_proteoforms;

			file_vertices_proteoforms.open(
					(path_disease_modules + trait + "/proteoformVertices.tsv").c_str(),
					std::ifstream::in);

			if (file_vertices_proteoforms.good()) {
				getline(file_vertices_proteoforms, proteoform);	// Skip the first line
				while (getline(file_vertices_proteoforms, proteoform, '\t')) {
					num_proteoforms++;

					string rest_of_line;
					getline(file_vertices_proteoforms, rest_of_line);

					std::smatch matches;
					std::regex expression { "[;,]\\d{5}" };
					std::regex_search(proteoform, matches, expression);

					if (matches.size() > 0) {
						num_modified++;
					}
				}
			} else {
				cout << "Failed to open the proteoforms file." << endl;
			}

			trait_to_num_modified[trait] = num_modified;
			ratio = num_proteoforms ?
					(double) num_modified / (double) num_proteoforms : 0.0;
			trait_to_ratio_modified[trait] = ratio;

			cout << trait << "\t" << num_proteoforms << "\t" << num_modified
					<< "\t" << ratio << endl;
			file_modified_percentage << trait << "\t" << num_proteoforms << "\t"
					<< num_modified << "\t" << ratio << endl;

			file_vertices_proteoforms.close();
		}
	}

	trait_file.close();
	file_modified_percentage.close();
}

set<string> getVertices(Entity entity, string trait) {

	// Read vertices file of the respective entity
	set<string> vertices;
	ifstream file_vertices;
	string value;				// Entity label: accession, gene name, etc.
	string line_leftover;

	file_vertices.open((path_disease_modules + trait + "/" + entity_str[entity] + "Vertices.tsv").c_str(), std::ifstream::in);

	if(file_vertices) {
		getline(file_vertices, line_leftover); 		// Discard the first line
		while(file_vertices >> value) {
			getline(file_vertices, line_leftover);
//			cout << " ++ " << value << endl;
			vertices.insert(value);
		}
	}
//	cout << "++++++++++" << endl;
	return vertices;
}

set<string> getOverlap(Entity entity, string one_trait, string other_trait) {

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

	create_ratios_file();

	// Read ratios file
	string trait;
	char trait_arr[80];
	int num_proteoforms = -1;
	int num_modified = -1;
	float ratio = -1.0;
	string path_file_modified_percentage = path_phenI + "modified_percentage.csv";
	multimap<string, string> candidates;
	ofstream file_candidates(path_phenI + "candidate_disease_pairs.csv");

	freopen(path_file_modified_percentage.c_str(), "r", stdin);

	// Read all traits
	while(scanf("%s%i%d%f", trait_arr, &num_proteoforms, &num_modified, &ratio) > 1) {
//		cout << trait_arr << " ** " << num_proteoforms << " ** " << num_modified << " ** " << ratio << endl;
		trait.assign(trait_arr);
		trait_to_ratio_modified[trait] = ratio;
	}

	fclose(stdin);

	file_candidates << "Trait1\t" << "Trait2\t" << "RatioModified1\t" << "RatioModified2\t" << "OverlapProteins\t" << "OverlapProteoforms" << endl;

	for(const auto &one_trait : trait_to_ratio_modified) {
		if (one_trait.second >= MIN_PERCENTAGE_MODIFIED) {											// For all phenotypes that have percentage higher than threshold
			for(const auto &other_trait : trait_to_ratio_modified) { 			// Compare to all other phenotypes for overlap
//				cout << "Comparing combination " << ++combination << ": " << one_trait.first << " with " << other_trait.first << endl;

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
