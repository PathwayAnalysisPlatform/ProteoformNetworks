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
#include <utility>

// Network (or graph with vertices and edges) containint all entities (genes, proteins, proteoforms and small molecules)
// in Reactome with all its interactions.

// Entities are sorted by category and assigned an index enumerating them alphabetically.
// Each entity type has a range of indexes [x, y], where all possible indexes between x and y inclusive are entities
// of the said type.
// First are the genes, then proteins, then proteoforms, then small molecules.
class Interactome {

    std::map<std::string, int> node_indexes;
    std::vector<std::string> node_names;
    std::vector<std::set<int>> adj_list;
//
//    std::vector<int> start_indexes;
//    std::vector<int> end_indexes;
//
//    void readRanges(std::string_view file_interactions);
//
//    void readEdges(std::string_view file_indexes);
//
//    std::map<int, std::vector<int>> genes_to_proteins;
//
//    std::map<int, std::vector<int>> proteins_to_proteoform;
//
//    void readGenesToProteins(std::string_view file_proteins_to_genes);
//
//    void readProteinsToProteoforms(std::string_view file_proteins_to_proteoforms);

public:

    Interactome();

    void addInteractions(std::vector<std::pair<int, int>> &interactions);

    std::vector<int> getNodes() const;

    [[nodiscard]] std::string getNodeName(int node) const;

//    Interactome(std::string_view file_vertices, std::string_view file_interactions, std::string_view file_ranges,
//                std::string_view file_proteins_to_genes, std::string_view file_proteins_to_proteoforms);
//
//    int index(std::string_view name);
//

//
//    Level get_type(int index);
//
//    Level get_type(std::string_view name);
//
//    bool isGene(int index);
//
//    bool isProtein(int index);
//
//    bool isProteoform(int index);
//
//    bool isSimpleEntity(int index);
//
//    std::vector<int> getProteins(int gene_index);
//
//    std::vector<int> getProteoforms(int protein_index);
//
//    std::set<int> getNeighbors(int index);
//
//    std::set<int> getSimpleEntityNeighbors(int index);
//
//    std::set<int> getProteinNeighbors(int index);
//
//    std::set<int> getProteoformNeighbors(int index);
//
//    std::vector<std::pair<int, int>> getInteractions(std::vector<int> indexes);
//
//    int getNumVertices();
//
//    int getStartIndexGenes();
//    int getEndIndexGenes();
//    int getStartIndexProteins();
//    int getEndIndexProteins();
//    int getStartIndexProteoforms();
//    int getEndIndexProteoforms();

    void addNode(int index);

    [[nodiscard]] std::vector<int> getInteractors(int node) const;

    [[nodiscard]] bool hasNode(int node) const;
    [[nodiscard]] bool hasNode(std::string_view name) const;

    void readNodeNames(std::istream &s);
};


#endif //PROTEOFORMNETWORKS_INTERACTOME_HPP
