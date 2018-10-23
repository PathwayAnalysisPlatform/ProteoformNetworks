#include "pathway_overlap.h"

const string path_file_proteoform_search = "resources/reactome/all_proteoforms/search.tsv";
const string file_pathway_pairs_with_key_ptms = "resources/reactome/pathway_pairs_with_key_ptms.txt";

int main(){
    set<pair<string, string>> pairs = findPairsWithKeyPTMExamples(0.0, path_file_proteoform_search, file_pathway_pairs_with_key_ptms);
}