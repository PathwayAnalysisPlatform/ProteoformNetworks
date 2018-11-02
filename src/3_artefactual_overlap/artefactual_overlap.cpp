#include <algorithm>
#include <bitset>
#include <ctime>
#include <exception>
#include <fstream>
#include <future>
#include <iostream>
#include <set>
#include <map>
#include <numeric>
#include <thread>
#include <unordered_set>
#include <vector>

using namespace std;

const size_t NUM_GENES = 23970;
const size_t NUM_PROTEINS = 10778;
const size_t NUM_PROTEOFORMS = 13911;

const int MIN_OVERLAP_SIZE = 2;
const int MAX_OVERLAP_SIZE = 10;

vector<string> loadEntities(const string& entities_file_path) {
    ifstream entities_file(entities_file_path);
    string entity, line_leftover;
    unordered_set<string> temp_set;
    vector<string> entities;

    if (!entities_file.is_open()) {
       throw std::runtime_error("Cannot open " + entities_file_path);
    }

    while (getline(entities_file, entity, '\t')) {
        temp_set.insert(entity);
        getline(entities_file, line_leftover);
    }
    entities.assign(temp_set.begin(), temp_set.end());
    sort(entities.begin(), entities.end());
    return entities;
}

map<string, int> fillMap(const vector<string>& index_to_entities) {
    map<string, int> entities_to_index;
    for (int I = 0; I < index_to_entities.size(); I++) {
        entities_to_index.emplace(index_to_entities[I], I);
    }
    return entities_to_index;
}

map<string, bitset<NUM_GENES>> loadPathwaysGeneMembers(const string& file_path, const map<string, int>& entities_to_index) {
    map<string, bitset<NUM_GENES>> result;
    ifstream file_search(file_path);
    string field, gene, pathway;
	bitset<NUM_GENES> empty_set;

    getline(file_search, field);                  // Skip csv header line
    while (getline(file_search, gene, '\t')) {  // Read the members of each pathway // Read GENE
        getline(file_search, field, '\t');        // Read UNIPROT
        getline(file_search, field, '\t');        // Read REACTION_STID
        getline(file_search, field, '\t');        // Read REACTION_DISPLAY_NAME
        getline(file_search, pathway, '\t');      // Read PATHWAY_STID
        getline(file_search, field);              // Read PATHWAY_DISPLAY_NAME

		if (result.find(pathway) == result.end()) {
			result.emplace(pathway, empty_set);
		}
		if (entities_to_index.find(gene) != entities_to_index.end()) {
         result[pathway].set(entities_to_index.find(gene)->second);
      }
    }
    return result;
}

map<string, bitset<NUM_PROTEINS>> loadPathwaysProteinMembers(const string& file_path, const map<string, int>& entities_to_index) {
	map<string, bitset<NUM_PROTEINS>> result;
	ifstream file_search(file_path);
	string field, entity, pathway;
	bitset<NUM_PROTEINS> empty_set;

	getline(file_search, field);                  // Skip csv header line
	while (getline(file_search, entity, '\t')) {  // Read the members of each pathway // Read UNIPROT
		getline(file_search, field, '\t');        // Read REACTION_STID
		getline(file_search, field, '\t');        // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');      // Read PATHWAY_STID
		getline(file_search, field);              // Read PATHWAY_DISPLAY_NAME

		if (result.find(pathway) == result.end()) {
			result.emplace(pathway, empty_set);
		}
		if (entities_to_index.find(entity) != entities_to_index.end()) {
         result[pathway].set(entities_to_index.find(entity)->second);
      }
	}
	return result;
}

map<string, bitset<NUM_PROTEOFORMS>> loadPathwaysProteoformMembers(const string& file_path, const map<string, int>& entities_to_index) {
	map<string, bitset<NUM_PROTEOFORMS>> result;
	ifstream file_search(file_path);
	string field, entity, pathway;
	bitset<NUM_PROTEOFORMS> empty_set;

	getline(file_search, field);                  // Skip csv header line
	while (getline(file_search, entity, '\t')) {  // Read the members of each pathway // Read PROTEOFORM
		getline(file_search, field, '\t');        // Read UNIPROT
		getline(file_search, field, '\t');        // Read REACTION_STID
		getline(file_search, field, '\t');        // Read REACTION_DISPLAY_NAME
		getline(file_search, pathway, '\t');      // Read PATHWAY_STID
		getline(file_search, field);              // Read PATHWAY_DISPLAY_NAME

		if (result.find(pathway) == result.end()) {
			result.emplace(pathway, empty_set);
		}
		if (entities_to_index.find(entity) != entities_to_index.end()) {
         result[pathway].set(entities_to_index.find(entity)->second);
      }
	}
	return result;
}

// SINGLE-CORE
template <size_t set_size>
set<pair<string, string>> findOverlappingPairs(const map<string, bitset<set_size>>& sets_to_members, int min_overlap, int max_overlap) {
	std::vector<typename map<string, bitset<set_size>>::const_iterator> nav;
   for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
      nav.push_back(it);
   }
std::cerr << "elementos: " << sets_to_members.size( ) << ", parejas: " << (sets_to_members.size( ) * sets_to_members.size( ) - sets_to_members.size( )) / 2 << "\n";
auto t0 = std::clock( );
   set<pair<string, string>> result;
	for (auto vit1 = nav.begin(); vit1 != nav.end(); vit1++) {
		for (auto vit2 = std::next(vit1); vit2 != nav.end(); vit2++) {
         bitset<set_size> overlap = (*vit1)->second & (*vit2)->second;
         if (min_overlap <= overlap.count() && overlap.count() <= max_overlap) {
            result.emplace((*vit1)->first, (*vit2)->first);
			}
		}
	}
auto t1 = std::clock( );
std::cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
	return result;
}

// MULTI-CORE
/*template <size_t set_size>
set<pair<string, string>> findOverlappingPairs(const map<string, bitset<set_size>>& sets_to_members, int min_overlap, int max_overlap) {
   std::vector<int> acc(sets_to_members.size( ));
   std::iota(acc.begin( ), acc.end( ), 0), std::partial_sum(acc.begin( ), acc.end( ), acc.begin( ));
   int pairs = (sets_to_members.size( ) * sets_to_members.size( ) - sets_to_members.size( )) / 2;
   int pairs_per_thread = pairs / std::thread::hardware_concurrency( ) + bool(pairs % std::thread::hardware_concurrency( ));
std::cerr << "elementos: " << sets_to_members.size( ) << ", parejas: " << pairs << "\n";
auto t0 = std::clock( );
	std::vector<typename map<string, bitset<set_size>>::const_iterator> nav;
   for (auto it = sets_to_members.begin(); it != sets_to_members.end(); it++) {
      nav.push_back(it);
   }
   std::vector<std::future<set<pair<string, string>>>> futures;
   for (int i = 0; i < std::thread::hardware_concurrency( ); ++i) {
      futures.push_back(std::async(std::launch::async, [&](int begin, int end) {
         set<pair<string, string>> result;
         for (int i = begin, row = 0, col; i < end; ++i) {
            while (acc[row] <= i) {
               ++row;
            }; col = i - acc[row - 1];
            bitset<set_size> overlap = nav[row]->second & nav[col]->second;
            if (min_overlap <= overlap.count() && overlap.count() <= max_overlap) {
               result.emplace(nav[row]->first, nav[col]->first);
            }
         }
         return result;
      }, i * pairs_per_thread, std::min(pairs, (i + 1) * pairs_per_thread)));
   }

   set<pair<string, string>> result;
   for (auto& fut : futures) {
      result.merge(fut.get( ));
   }
auto t1 = std::clock( );
std::cerr << "tardamos " << double(t1 - t0) / CLOCKS_PER_SEC << "\n";
	return result;
}*/

template <size_t set_size>
set<pair<string, string>> findOverlappingPairs(const map<string, bitset<set_size>>& sets_to_members) {
    return findOverlappingPairs(sets_to_members, 1, static_cast<int>(set_size));
}

/**
 * Find artefactual overlaps: pairs of pathways that share nodes only in the
 * gene or protein level, but not at the proteoform level
 */
set<pair<string, string>> findPathwayPairsWithArtefactualOverlap(const map<string, bitset<NUM_GENES>>& pathway_to_genes,
                                                                 const map<string, bitset<NUM_PROTEINS>>& pathway_to_proteins,
                                                                 const map<string, bitset<NUM_PROTEOFORMS>>& pathway_to_proteoforms,
                                                                 const string& path_file_gene_art_pairs,
                                                                 const string& path_file_protein_art_pairs) {
    set<pair<string, string>> result;
    set<pair<string, string>> overlaping_gene_set_pairs;
    set<pair<string, string>> overlaping_protein_set_pairs;
    set<pair<string, string>> overlaping_proteoform_set_pairs;
    ofstream file_gene_art_pairs(path_file_gene_art_pairs);
    ofstream file_protein_art_pairs(path_file_protein_art_pairs);

    cout << "Calculating gene sets overlap..." << endl;
    overlaping_gene_set_pairs = findOverlappingPairs(pathway_to_genes, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE);

    cout << "Calculating protein sets overlap..." << endl;
    overlaping_protein_set_pairs = findOverlappingPairs(pathway_to_proteins, MIN_OVERLAP_SIZE, MAX_OVERLAP_SIZE);

    cout << "Calculating proteoform sets..." << endl;
    overlaping_proteoform_set_pairs = findOverlappingPairs(pathway_to_proteoforms);

    cout << "Comparing gene and proteoform pairs..." << endl;
    for (const auto &gene_set_pair : overlaping_gene_set_pairs) {
        set<pair<string, string>>::iterator it = overlaping_proteoform_set_pairs.find(gene_set_pair);
        if (it != overlaping_proteoform_set_pairs.end()) {
            cout << gene_set_pair.first << "\t" << gene_set_pair.second << "\n";
            file_gene_art_pairs << gene_set_pair.first << "\t" << gene_set_pair.second << "\n";
            result.insert(gene_set_pair);
        }
    }

    cout << "Comparing protein and proteoform pairs..." << endl;
    for (const auto &protein_set_pair : overlaping_protein_set_pairs) {
        set<pair<string, string>>::iterator it = overlaping_proteoform_set_pairs.find(protein_set_pair);
        if (it != overlaping_proteoform_set_pairs.end()) {
            cout << protein_set_pair.first << "\t" << protein_set_pair.second << "\n";
            file_protein_art_pairs << protein_set_pair.first << "\t" << protein_set_pair.second << "\n";
            result.insert(protein_set_pair);
        }
    }

    return result;
}

// Find pairs of modules/pathways that overlap on gene or protein network, but not in proteoform network
void doArtefactualOverlapAnalysis(const string& path_file_gene_search,
	const string& path_file_protein_search,
	const string& path_file_proteoform_search,
	const string& report_path,
	const string& path_file_gene_art_pairs,
	const string& path_file_protein_art_pairs) {
	ofstream report(report_path);

	// TODO: Report pairs of disease modules
	//findDiseaseModulePairsWithArtefactualOverlap(path_file_gene_art_pairs, path_file_protein_art_pairs);

	vector<string> index_to_genes = loadEntities(path_file_gene_search);
	map<string, int> genes_to_index = fillMap(index_to_genes);

	vector<string> index_to_proteins = loadEntities(path_file_protein_search);
	map<string, int> proteins_to_index = fillMap(index_to_proteins);

	vector<string> index_to_proteoforms = loadEntities(path_file_proteoform_search);
	map<string, int> proteoforms_to_index = fillMap(index_to_proteoforms);

	map<string, bitset<NUM_GENES>> pathway_to_genes = loadPathwaysGeneMembers(path_file_gene_search, genes_to_index);
	map<string, bitset<NUM_PROTEINS>> pathway_to_proteins = loadPathwaysProteinMembers(path_file_protein_search, proteins_to_index);
	map<string, bitset<NUM_PROTEOFORMS>> pathway_to_proteoforms = loadPathwaysProteoformMembers(path_file_proteoform_search, proteoforms_to_index);

	set<pair<string, string>> examples = findPathwayPairsWithArtefactualOverlap(pathway_to_genes,
		pathway_to_proteins,
		pathway_to_proteoforms,
		path_file_gene_art_pairs,
		path_file_protein_art_pairs);

	report << "PATHWAY_1\tPATHWAY_2\tOVERLAP_GENES\tOVERLAP_PROTEINS\tOVERLAP_PROTEOFORMS";
	for (const auto &example : examples) {
		report << example.first << " " << example.second << "\t";

		bitset<NUM_GENES> overlap_genes = pathway_to_genes[example.first] & pathway_to_genes[example.second];
		int printed = 0;
		for (int I = 0; I < NUM_GENES; I++) {
			if (overlap_genes.test(I)) {
				report << index_to_genes[I];
				printed++;
				if (printed != overlap_genes.count()) {
					report << ";";
				}
			}
		}
		report << "\t";

		bitset<NUM_PROTEINS> overlap_proteins = pathway_to_proteins[example.first] & pathway_to_proteins[example.second];
		printed = 0;
		for (int I = 0; I < NUM_PROTEINS; I++) {
			if (overlap_proteins.test(I)) {
				report << index_to_proteins[I];
				printed++;
				if (printed != overlap_proteins.count()) {
					report << ";";
				}
			}
		}
		report << "\t";
		report << "\n";
	}
}

/* set<pair<string, string>> findDiseaseModulePairsWithArtefactualOverlap(string path_file_gene_art_pairs,
                                                                       string path_file_protein_art_pairs) {
    CreatePhenotypePairsFile(0.95, protein);

    multimap<string, string> pairs;
    ifstream pairs_file("../resources/PheGenI/pairs_diff_overlap_genes_vs_proteoforms.csv");

    if (pairs_file) {
        string one_trait, other_trait;
        string line_leftover;
        while (pairs_file >> one_trait >> other_trait) {
            getline(pairs_file, line_leftover);
            pairs.insert(make_pair(one_trait, other_trait));
        }
    }

    pairs.insert(make_pair("Aorta", "Asthma"));
    pairs.insert(make_pair("ColitisUlcerative", "Creatinine"));
    pairs.insert(make_pair("Creatinine", "Adiponectin"));
    pairs.insert(make_pair("Creatinine", "Prostate_SpecificAntigen"));
    pairs.insert(make_pair("Creatinine", "WetMacularDegeneration"));
    pairs.insert(make_pair("MacularDegeneration", "Hypothyroidism"));
    pairs.insert(make_pair("WaistCircumference", "Hypothyroidism"));
    pairs.insert(make_pair("WetMacularDegeneration", "Hypothyroidism"));
    ShowOverlapsAndSplittingProteoforms(pairs);
} */
