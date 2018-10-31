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


const int MIN_OVERLAP_SIZE = 0;
const int MAX_OVERLAP_SIZE = 2000;

const std::regex RGX_ACCESSION_DELIMITER{ "[;-]" };
const std::regex RGX_MODIFICATION{ "[;,]\\d{5}" };

bool isModified(const string& proteoform) {
	std::smatch modification;
	return std::regex_search(proteoform, modification, RGX_MODIFICATION);
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
set<pair<string, string>> findPairsWithKeyPTMExamples(float min_modified_percentage, string path_file_proteoform_search,
	string path_file_pathways_in_key_ptm_overlaps, string path_file_proteins_in_key_ptm_overlaps,
	string path_file_proteoforms_in_key_ptm_overlaps, string path_file_modifications_in_key_ptm_overlaps)
{
	ifstream file_proteoform_search(path_file_proteoform_search);

	ofstream file_pathway_pairs_in_key_ptm_overlaps(path_file_pathways_in_key_ptm_overlaps);
	ofstream file_proteins_in_key_ptm_overlaps(path_file_proteins_in_key_ptm_overlaps);
	ofstream file_proteoforms_in_key_ptm_overlaps(path_file_proteoforms_in_key_ptm_overlaps);
	ofstream file_modifications_in_key_ptm_overlaps(path_file_modifications_in_key_ptm_overlaps);

	string field, pathway, proteoform;
	map<string, set<string>> pathways_to_proteoforms;
	map<string, set<string>> proteoforms_to_pathways;
	multimap<string, string> candidate_pairs;
	set<string> proteoforms;
	map<string, int> proteins_freq;		// Times a protein is in the overlap of a pathway pair
	map<string, int> modifications_freq;	// Times a modification is in the overlap of a pathway pair

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

		proteoforms.insert(proteoform);
		proteoforms_to_pathways[proteoform].insert(pathway);
		pathways_to_proteoforms[pathway].insert(proteoform);
	}

	file_proteins_in_key_ptm_overlaps << "PROTEIN\tFREQUENCY\n";
	file_proteoforms_in_key_ptm_overlaps << "PROTEOFORM\tFREQUENCY\n";
	file_modifications_in_key_ptm_overlaps << "MODIFICATION\tFREQUENCY\n";
	file_pathway_pairs_in_key_ptm_overlaps << "PATHWAY1\tPATHWAY2\tOVERLAP_SIZE\tMODIFIED_RATIO\n";

	// Check which modified proteoforms appear in more than one pathway
	for (const auto &proteoform : proteoforms) {
		int num_pathways = proteoforms_to_pathways[proteoform].size();
		if (isModified(proteoform) && num_pathways > 1) {
			int num_pathway_pairs = ((num_pathways * num_pathways) - num_pathways) / 2;
			cout << "Found candidate proteoform: " << proteoform << " in " << num_pathways << " pathways and " << num_pathway_pairs << " pairs" << "\n";
																				// Print proteoform frequencies
			file_proteoforms_in_key_ptm_overlaps << proteoform << "\t" << num_pathway_pairs << "\n";

			proteins_freq[getAccession(proteoform)] += num_pathway_pairs;			// Update protein frequencies
			for (const auto &modification : getModifications(proteoform)) {		// Update modification frequencies
				modifications_freq[modification] += num_pathway_pairs;
			}

			auto ordered_pairs = calculateOrderedPairs(proteoforms_to_pathways[proteoform]);
			candidate_pairs.insert(ordered_pairs.begin(), ordered_pairs.end());
		}
	}

	for (const auto &pair : proteins_freq) {									// Print protein frequencies
		file_proteins_in_key_ptm_overlaps << pair.first << "\t" << pair.second << "\n";
	}

	for (const auto &pair : modifications_freq) {								// Print modification frequencies
		file_modifications_in_key_ptm_overlaps << pair.first << "\t" << pair.second << "\n";
	}

	// Check all pairs of candidate pathways
	for (const auto &candidate_pair : candidate_pairs) {
		set<string> overlap = calculateOverlap(candidate_pair.first, candidate_pair.second, pathways_to_proteoforms);
		if (overlap.size() > 0) {
			float ratio = calculateModifiedRatio(overlap);
			if (ratio > min_modified_percentage)
			{
				cout << "Found example: " << candidate_pair.first << " with " << candidate_pair.second << "\n";
				result.emplace(candidate_pair.first, candidate_pair.second);
				file_pathway_pairs_in_key_ptm_overlaps 
					<< candidate_pair.first << "\t" 
					<< candidate_pair.second << "\t" 
					<< overlap.size() << "\t" 
					<< ratio << "\n";

				/*file_pathway_pairs_in_key_ptm_overlaps << "===============================================\n";
				file_pathway_pairs_in_key_ptm_overlaps << candidate_pair.first << " with " << candidate_pair.second << "\n";
				file_pathway_pairs_in_key_ptm_overlaps << "-----------------------------------------------\n";
				for (const auto &member : overlap) {
					file_pathway_pairs_in_key_ptm_overlaps << member << "\t";
				}
				file_pathway_pairs_in_key_ptm_overlaps << "\n";
				file_pathway_pairs_in_key_ptm_overlaps << "===============================================\n";*/
			}
		}
	}
	return result;
}

string getAccession(string proteoform) {
	smatch match_end_of_accession;
	if (!regex_search(proteoform, match_end_of_accession, RGX_ACCESSION_DELIMITER)) {
		return proteoform;
	}
	return proteoform.substr(0, match_end_of_accession.position(0));
}

set<string> getModifications(string proteoform) {
	set<string> modifications;
	sregex_iterator it(proteoform.begin(), proteoform.end(), RGX_MODIFICATION);
	sregex_iterator end;
	while (it != end) {
		if (it->str().find(';') || it->str().find(',')) {
			modifications.insert(it->str().substr(1));
		}
		else {
			modifications.insert((*it)[0]);
		}
		it++;
	}
	return modifications;
}

set<pair<string, string>> calculateOrderedPairs(set<string> set_of_ids) {
	set<pair<string, string>> result;
	for (const auto &one_id : set_of_ids) {
		for (const auto &other_id : set_of_ids) {
			if (one_id.compare(other_id) < 0) {
				result.emplace(one_id, other_id);
			}
		}
	}

	return result;
}

/**
 * Find pairs of pathways that only share nodes
 */
set<pair<string, string>> findOverlappingPairs(string path_file_search, int min_overlap_size, int max_overlap_size) {

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

			if (overlap.size() >= min_overlap_size && overlap.size() <= max_overlap_size) {
				result.emplace(one_pathway, other_pathway);
				// cout << ".";
			}
		}
	}

	return result;
}

set<pair<string, string>> findOverlappingPairs(string path_file_search) {
	return findOverlappingPairs(path_file_search, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE);
}

map < string, set<string>> calculateMembers(string path_file_search) {
	ifstream file_search(path_file_search);
	string field, pathway, entity;
	map<string, set<string>> pathways_to_proteoforms;

	if (!file_search.good()) {
		cout << "Could not open search file.\n";
	}

	getline(file_search, field); // Skip csv header line
	while (getline(file_search, entity, '\t'))
	{													// Read the members of each pathway
		getline(file_search, field, '\t');   // Read UNIPROT
		getline(file_search, field, '\t');   // Read REACTION_STID
		getline(file_search, field, '\t');   // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t'); // Read PATHWAY_STID
		getline(file_search, field);			// Read PATHWAY_DISPLAY_NAME

		pathways_to_proteoforms[pathway].insert(entity);
	}
	return pathways_to_proteoforms;
}

set<string> calculateOverlap(string one_pathway, string other_pathway, map<string, set<string>> pathway_to_members) {
	set<string> overlap;
	set_intersection(
		pathway_to_members[one_pathway].begin(), pathway_to_members[one_pathway].end(),
		pathway_to_members[other_pathway].begin(), pathway_to_members[other_pathway].end(),
		inserter(overlap, overlap.begin()));
	return overlap;
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
	map<string, set<string>> pathways_to_proteoforms;

	cout << "Calculating gene sets overlap..." << endl;
	sets_genes = findOverlappingPairs(path_file_gene_search, 2, 10);

	cout << "Calculating protein sets overlap..." << endl;
	sets_proteins = findOverlappingPairs(path_file_protein_search, 2, 10);

	cout << "Calculating proteoform sets..." << endl;
	pathways_to_proteoforms = calculateMembers(path_file_proteoform_search);

	cout << "Comparing gene and proteoform pairs..." << endl;
	ofstream file_gene_art_pairs(path_file_gene_art_pairs);
	for (const auto &pair : sets_genes) {
		if (calculateOverlap(pair.first, pair.second, pathways_to_proteoforms).size() == 0) {
			cout << pair.first << "\t" << pair.second << "\n";
			file_gene_art_pairs << pair.first << "\t" << pair.second << "\n";
			result.emplace(pair.first, pair.second);
		}
	}

	cout << "Comparing protein and proteoform pairs..." << endl;
	ofstream file_protein_art_pairs(path_file_protein_art_pairs);
	for (const auto &pair : sets_genes) {
		if (calculateOverlap(pair.first, pair.second, pathways_to_proteoforms).size() == 0) {
			cout << pair.first << "\t" << pair.second << "\n";
			file_gene_art_pairs << pair.first << "\t" << pair.second << "\n";
			result.emplace(pair.first, pair.second);
		}
	}

	return result;
}