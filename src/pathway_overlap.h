#ifndef PATHWAY_OVERLAP_H_
#define PATHWAY_OVERLAP_H_

#include <set>
#include <map>
#include <string>
#include <cstdio>

using namespace std;

set<pair<string, string>> findPathwayPairsWithKeyPTMExamples(
		float min_modified_percentage, string path_file_proteoform_search);

set<pair<string, string>> findPathwayPairsWithArtifactualOverlapExamples(
		double minModifiedPercentage, string path_gene_search_file,
		string path_protein_search_file, string path_file_proteoform_search);

#endif /* PATHWAY_OVERLAP_H_ */
