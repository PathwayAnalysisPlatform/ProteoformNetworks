//============================================================================
// Description : Get set pairs with different overlap for genes, proteins
//               and proteoforms.
// Author      : Luis Francisco Hern�ndez S�nchez
// Copyright   : Licence Apache 2.0
//============================================================================

#include "phenotype_overlap.h"
#include "pathway_overlap.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

using namespace std;

// Input files
const string path_file_gene_search = "resources/reactome/all_genes/search.tsv";
const string path_file_protein_search = "resources/reactome/all_proteins/search.tsv";
string path_file_proteoform_search = "resources/reactome/all_proteoforms/test.txt";

// Output files
const string path_file_pathways_in_key_ptm_overlaps = "resources/reactome/key_ptm_overlap/pathway_pairs.txt";
const string path_file_proteins_in_key_ptm_overlaps = "resources/reactome/key_ptm_overlap/proteins.txt";
const string path_file_proteoforms_in_key_ptm_overlaps = "resources/reactome/key_ptm_overlap/proteoforms.txt";
const string path_file_modifications_in_key_ptm_overlaps = "resources/reactome/key_ptm_overlap/modifications.txt";

int main(int argc, char *argv[]) {
	if (argc >= 1) {
		path_file_proteoform_search = argv[1];
	}

	set<pair<string, string>> pairs;

	pairs = findPairsWithKeyPTMExamples( 
		0.0, 
		path_file_proteoform_search,
		path_file_pathways_in_key_ptm_overlaps,
		path_file_proteins_in_key_ptm_overlaps,
		path_file_proteoforms_in_key_ptm_overlaps,
		path_file_modifications_in_key_ptm_overlaps);


	/*pairs = findPathwayPairsWithArtifactualOverlapExamples(
		path_file_gene_search, 
		path_file_protein_search, 
		path_file_proteoform_search, 
		"resources/reactome/gene_art_pairs.txt", 
		"resources/reactome/protein_art_pairs.txt");*/

	//CreatePhenotypePairsFile(0.95, protein);
	// Read all pairs with opo and sep types

	/*multimap<string, string> pairs;
	 ifstream pairs_file("../resources/PheGenI/pairs_diff_overlap_genes_vs_proteoforms.csv");

	 if(pairs_file) {
	 string one_trait, other_trait;
	 string line_leftover;
	 while(pairs_file >> one_trait >> other_trait) {
	 getline(pairs_file, line_leftover);
	 pairs.insert(make_pair(one_trait, other_trait));
	 }
	 }*/

//	pairs.insert(make_pair("Aorta", "Asthma"));
//	pairs.insert(make_pair("ColitisUlcerative", "Creatinine"));
//	pairs.insert(make_pair("Creatinine", "Adiponectin"));
//	pairs.insert(make_pair("Creatinine", "Prostate_SpecificAntigen"));
//	pairs.insert(make_pair("Creatinine", "WetMacularDegeneration"));
//	pairs.insert(make_pair("MacularDegeneration", "Hypothyroidism"));
//	pairs.insert(make_pair("WaistCircumference", "Hypothyroidism"));
//	pairs.insert(make_pair("WetMacularDegeneration", "Hypothyroidism"));
	//ShowOverlapsAndSplittingProteoforms(pairs);
	return 0;
}
