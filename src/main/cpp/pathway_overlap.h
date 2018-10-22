#ifndef PATHWAY_OVERLAP_H_
#define PATHWAY_OVERLAP_H_

#include <set>
#include <map>
#include <string>
#include <cstdio>

using namespace std;

set<pair<string, string>> findPairsWithKeyPTMExamples(float min_modified_percentage, string path_file_proteoform_search, string path_file_pathway_pairs_with_key_ptms);

set<pair<string, string>> findPathwayPairsWithArtifactualOverlapExamples(string path_file_gene_search, string path_file_protein_search, string path_file_proteoform_search,
	string path_file_gene_art_pairs, string path_file_protein_art_pairs);


#endif /* PATHWAY_OVERLAP_H_ */
