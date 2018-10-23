#ifndef PATHWAY_OVERLAP_H_
#define PATHWAY_OVERLAP_H_

#include <set>
#include <map>
#include <string>
#include <cstdio>

using namespace std;

set<pair<string, string>> findPairsWithKeyPTMExamples(float min_modified_percentage, string path_file_proteoform_search,
	string path_file_pathway_pairs_in_key_ptm_overlaps, string path_file_proteins_in_key_ptm_overlaps,
	string path_file_proteoforms_in_key_ptm_overlaps, string path_file_modifications_in_key_ptm_overlaps);

set<pair<string, string>> findPathwayPairsWithArtifactualOverlapExamples(
	string path_file_gene_search, string path_file_protein_search, string path_file_proteoform_search,
	string path_file_gene_art_pairs, string path_file_protein_art_pairs);

set<string> calculateOverlap(string one_pathway, string other_pathway, map<string, set<string>> pathway_to_members);

string getAccession(string proteoform);

set<string> getModifications(string proteoform);

set<pair<string, string>> calculateOrderedPairs(set<string> set_of_ids);

#endif /* PATHWAY_OVERLAP_H_ */