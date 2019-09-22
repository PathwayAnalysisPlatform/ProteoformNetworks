#include "reactome.hpp"

using namespace std;

umsb loadGeneSets(string_view file_path, const bimap_str_int &phegeni_genes, bool pathways) {
    umsb result;
    // TODO: Initialize result structure to the bitsets
    std::ifstream file_search(file_path.data());
    std::string field, gene, reaction, pathway;

    getline(file_search, field);                // Skip csv header line
    while (getline(file_search, gene, '\t')) {  // Read the members of each gene_set // Read GENE
        getline(file_search, field, '\t');       // Read UNIPROT
        getline(file_search, reaction, '\t');    // Read REACTION_STID
        getline(file_search, field, '\t');       // Read REACTION_DISPLAY_NAME
        getline(file_search, pathway, '\t');     // Read PATHWAY_STID
        getline(file_search, field);             // Read PATHWAY_DISPLAY_NAME

        result[pathways ? pathway : reaction][phegeni_genes.str_to_int.find(gene)->second].set();
    }
    return result;
}


umss loadPathwayNames(string_view path_search_file) {
    umss result;
    ifstream file_search(path_search_file.data());
    string field, pathway, pathway_name;
    int field_cont = 0;
    int pathway_stid_column = 0;

    while (getline(file_search, field, '\t')) {
        if (field == "PATHWAY_STID") {
            pathway_stid_column = field_cont;
            break;
        }
        field_cont++;
    }
    getline(file_search, field);  // Skip header line leftoever

    while (getline(file_search, field, '\t')) {                        // While there are records in the file
        for (int column = 1; column < pathway_stid_column; column++) {  // Skip fields before the pathway_stid
            getline(file_search, field, '\t');
        }
        getline(file_search, pathway, '\t');  // Read PATHWAY_STID
        getline(file_search, pathway_name);   // Read PATHWAY_DISPLAY_NAME
        result[pathway] = pathway_name;
    }
    return result;
}