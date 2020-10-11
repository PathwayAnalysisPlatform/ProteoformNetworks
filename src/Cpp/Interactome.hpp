#ifndef PROTEOFORMNETWORKS_INTERACTOME_HPP
#define PROTEOFORMNETWORKS_INTERACTOME_HPP

#include <string_view>
#include <map>
#include <vector>
#include <set>
#include "bimap_str_int.hpp"
#include "types.hpp"
#include <fstream>
#include <iostream>

class Interactome {
    Bimap_str_int vertices;
    std::map<int, std::set<int>> adj;

    std::vector<int> start_indexes;
    std::vector<int> end_indexes;

    void readInteractions(std::string_view file_interactions);

    void readIndexes(std::string_view file_indexes);

    std::map<int, std::vector<int>> genes_to_proteins;

    std::map<int, std::vector<int>> proteins_to_proteoform;

    void readGenesToProteins(std::string_view file_proteins_to_genes);

    void readProteinsToProteoforms(std::string_view file_proteins_to_proteoforms);

public:

    Interactome(std::string_view file_vertices, std::string_view file_interactions, std::string_view file_indexes,
                std::string_view file_proteins_to_genes, std::string_view file_proteins_to_proteoforms);

    int index(std::string_view name);

    Level get_type(int index);

    Level get_type(std::string_view name);

    bool isGene(int index);

    bool isProtein(int index);

    bool isProteoform(int index);

    bool isSimpleEntity(int index);

    std::vector<int> getProteins(int gene_index);

    std::vector<int> getProteoforms(int protein_index);

    std::set<int> getNeighbors(int index);

    std::set<int> getSimpleEntityNeighbors(int index);

    std::set<int> getProteinNeighbors(int index);

    std::set<int> getProteoformNeighbors(int index);

    std::vector<std::pair<int, int>> getInteractions(std::vector<int> indexes);

};


#endif //PROTEOFORMNETWORKS_INTERACTOME_HPP
