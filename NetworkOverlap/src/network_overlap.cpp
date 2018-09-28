//============================================================================
// Description : Get pathway pairs with different overlap for genes, proteins
//               and proteoforms.
// Author      : Luis Francisco Hernández Sánchez
// Copyright   : Licence Apache 2.0
//============================================================================

#include "phenotype_overlap.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

using namespace std;

int main() {

	CreatePhenotypePairsFile(0.25, protein);
	// Read all pairs with opo and sep types

	multimap<string, string> pairs;
	ifstream pairs_file("../resources/PheGenI/pairs_diff_overlap_genes_vs_proteoforms.csv");

	if(pairs_file) {
		string one_trait, other_trait;
		string line_leftover;
		while(pairs_file >> one_trait >> other_trait) {
			getline(pairs_file, line_leftover);
			pairs.insert(make_pair(one_trait, other_trait));
		}
	}

//	pairs.insert(make_pair("Aorta", "Asthma"));
//	pairs.insert(make_pair("ColitisUlcerative", "Creatinine"));
//	pairs.insert(make_pair("Creatinine", "Adiponectin"));
//	pairs.insert(make_pair("Creatinine", "Prostate_SpecificAntigen"));
//	pairs.insert(make_pair("Creatinine", "WetMacularDegeneration"));
//	pairs.insert(make_pair("MacularDegeneration", "Hypothyroidism"));
//	pairs.insert(make_pair("WaistCircumference", "Hypothyroidism"));
//	pairs.insert(make_pair("WetMacularDegeneration", "Hypothyroidism"));

	ShowOverlapsAndSplittingProteoforms(pairs);

	return 0;
}
