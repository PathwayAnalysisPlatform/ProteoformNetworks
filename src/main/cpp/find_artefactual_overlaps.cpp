#include "pathway_overlap.h"

const string path_file_gene_search = "resources/reactome/all_genes/search.tsv";
const string path_file_protein_search = "resources/reactome/all_proteins/search.tsv";
const string path_file_proteoform_search = "resources/reactome/all_proteoforms/search.tsv";

const string path_file_gene_art_pairs = "resources/reactome/gene_art_pairs.txt";
const string path_file_protein_art_pairs = "resources/reactome/protein_art_pairs.txt";

int main() {
findPathwayPairsWithArtifactualOverlapExamples(
		path_file_gene_search, 
		path_file_protein_search, 
		path_file_proteoform_search, 
		path_file_gene_art_pairs, 
		path_file_protein_art_pairs);
}