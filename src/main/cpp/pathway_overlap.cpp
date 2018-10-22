//============================================================================
// Description : Get pathway pairs that overlap only with modified proteins,
// 				 or that artefactually overlap at the gene or protein level
// Author      : Luis Francisco Hern�ndez S�nchez
// Copyright   : Licence Apache 2.0
//============================================================================

#include "pathway_overlap.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <string>

using namespace std;

const std::regex expression{ "[;,]\\d{5}" };

bool isModified(const string& proteoform) {
	std::smatch modification;
	return std::regex_search(proteoform, modification, expression);
}

bool areAllModified(set<string> proteoforms) {
	for (const auto &proteoform : proteoforms) {
		if (!isModified(proteoform)) {
			return false;
		}
	}
	return true;
}

float calculateModifiedRatio(set<string> proteoforms) {
	int cont = 0;

	for (const auto &proteoform : proteoforms) {
		if (isModified(proteoform)) {
			cont++;
		}
	}

	return (float)cont / (float)proteoforms.size();
}

/**
 * Find key PTM examples: pairs of pathways or reactions that only share
 * members that are proteins with at least one modification
 */
set<pair<string, string>> findPairsWithKeyPTMExamples(float min_modified_percentage, string path_file_proteoform_search, string path_file_pathway_pairs_with_key_ptms)
{
	ifstream file_proteoform_search(path_file_proteoform_search);
	ofstream file_pathway_pairs_with_key_ptms(path_file_pathway_pairs_with_key_ptms);
	string field, pathway, proteoform;
	set<string> proteoforms;
	map<string, set<string>> pathways_to_proteoforms;
	map<string, set<string>> proteoforms_to_pathways;
	multimap<string, string> candidate_pairs;
	set<pair<string, string>> result; // Insert always pairs ordering lexicographically the members

	if (!file_proteoform_search.good()) {
		cout << "Could not open search file.\n";
	}

	getline(file_proteoform_search, field); // Skip csv header line
	while (getline(file_proteoform_search, proteoform, '\t'))
	{													// Read the members of each pathway
		getline(file_proteoform_search, field, '\t');   // Read UNIPROT
		getline(file_proteoform_search, field, '\t');   // Read REACTION_STID
		getline(file_proteoform_search, field, '\t');   // Read REACTION_DISPLAY_NAME
		getline(file_proteoform_search, pathway, '\t'); // Read PATHWAY_STID
		getline(file_proteoform_search, field);			// Read PATHWAY_DISPLAY_NAME

		if (isModified(proteoform))
		{
			proteoforms.insert(proteoform);
		}
		proteoforms_to_pathways[proteoform].insert(pathway);
		pathways_to_proteoforms[pathway].insert(proteoform);
	}

	cout << "Finished finding candidate proteoforms.\n";

	// Check which modified proteoforms appear in more than one pathway

	for (const auto &proteoform : proteoforms) {
		// cout << "Checking " << proteoform << "\n";
		set<string> containing_pathways = proteoforms_to_pathways[proteoform];
		if (containing_pathways.size() > 0) {
			// cout << "Found candidate proteoform: " << proteoform << "\n";
			for (const auto &one_pathway : containing_pathways) {
				for (const auto &other_pathway : containing_pathways) {
					if (one_pathway.compare(other_pathway) < 0) {
						candidate_pairs.emplace(one_pathway, other_pathway);
					}
				}
			}
		}
	}

	// Check all pairs of candidate pathways
	for (const auto &candidate_pair : candidate_pairs)
	{
		
		set<string> overlap;
		set<string> one_set = pathways_to_proteoforms[candidate_pair.first];
		set<string> other_set = pathways_to_proteoforms[candidate_pair.second];
		set_intersection(one_set.begin(), one_set.end(),
			other_set.begin(), other_set.end(),
			inserter(overlap, overlap.begin()));

		if (overlap.size() > 0)
		{
			if (calculateModifiedRatio(overlap) > min_modified_percentage)
			{
				cout << "Found key ptms at " << candidate_pair.first << " with " << candidate_pair.second << "\n";
				result.emplace(candidate_pair.first, candidate_pair.second);
				file_pathway_pairs_with_key_ptms << "===============================================\n";
				file_pathway_pairs_with_key_ptms << candidate_pair.first << " with " << candidate_pair.second << "\n";
				file_pathway_pairs_with_key_ptms << "-----------------------------------------------\n";
				for (const auto &member : overlap)
				{
					file_pathway_pairs_with_key_ptms << member << "\t";
				}
				file_pathway_pairs_with_key_ptms << "\n";
				file_pathway_pairs_with_key_ptms << "===============================================\n";
			}
		}
	}
	return result;
}

/**
 * Find pairs of pathways that only share nodes
 */
set<pair<string, string>> findOverlappingPairs(string path_file_search) {

	ifstream file_search(path_file_search);
	string field, pathway, entity;
	set<string> pathways;
	map<string, set<string>> pathways_to_entities;
	set<pair<string, string>> result;					// Insert always pairs ordering lexicographically the members

	getline(file_search, field);						// Skip csv header line
	while (getline(file_search, entity, '\t')) {		// Read the members of each pathway
		getline(file_search, field, '\t'); 				// Read UNIPROT
		getline(file_search, field, '\t');				// Read REACTION_STID
		getline(file_search, field, '\t');				// Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');			// Read PATHWAY_STID
		getline(file_search, field);					// Read PATHWAY_DISPLAY_NAME

		pathways.insert(pathway);
		pathways_to_entities[pathway].insert(entity);
	}

	for (const auto &one_pathway : pathways) {
		set<string> one_set = pathways_to_entities[one_pathway];

		for (const auto &other_pathway : pathways) {
			if (one_pathway.compare(other_pathway) >= 0) {
				continue;
			}

			set<string> other_set = pathways_to_entities[other_pathway];
			set<string> overlap;
			set_intersection(one_set.begin(), one_set.end(), other_set.begin(),
				other_set.end(), inserter(overlap, overlap.begin()));

			if (overlap.size() > 0) {
				result.emplace(one_pathway, other_pathway);
				// cout << ".";
			}
		}
	}

	return result;
}

/**
 * Find artefactual overlaps: pairs of pathways that share nodes only in the
 * gene or protein level, but not at the proteoform level
 */
set<pair<string, string>> findPathwayPairsWithArtifactualOverlapExamples(string path_file_gene_search, string path_file_protein_search, string path_file_proteoform_search,
	string path_file_gene_art_pairs, string path_file_protein_art_pairs) {
	set<pair<string, string>> result;			// Insert always pairs ordering lexicographically the members

	set<pair<string, string>> sets_genes;
	set<pair<string, string>> sets_proteins;
	set<pair<string, string>> sets_proteoforms;

	cout << "Calculating gene sets overlap..." << endl;
	sets_genes = findOverlappingPairs(path_file_gene_search);

	cout << "Calculating protein sets overlap..." << endl;
	sets_proteins = findOverlappingPairs(path_file_protein_search);

	cout << "Calculating proteoform sets overlap..." << endl;
	sets_proteoforms = findOverlappingPairs(path_file_proteoform_search);

	cout << "Comparing gene and proteoform pairs..." << endl;
	ofstream file_gene_art_pairs(path_file_gene_art_pairs);
	for (const auto &pair : sets_genes) {
		if (sets_proteoforms.find(pair) != sets_proteoforms.end()) {
			cout << pair.first << "\t" << pair.second << "\n";
			file_gene_art_pairs << pair.first << "\t" << pair.second << "\n";
		}
	}

	cout << "Comparing protein and proteoform pairs..." << endl;
	ofstream file_protein_art_pairs(path_file_protein_art_pairs);
	for (const auto &pair : sets_genes) {
		if (sets_proteoforms.find(pair) != sets_proteoforms.end()) {
			cout << pair.first << "\t" << pair.second << "\n";
			file_protein_art_pairs << pair.first << "\t" << pair.second << "\n";
		}
	}

	return result;
}