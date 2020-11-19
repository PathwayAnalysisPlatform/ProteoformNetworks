#ifndef PROTEOFORMNETWORKS_MODULE_HPP
#define PROTEOFORMNETWORKS_MODULE_HPP


#include <string>
#include <map>
#include <set>
#include <vector>
#include "Interactome.hpp"
#include "types.hpp"

class Module {
    std::string name;
    Level level;
    std::map<int, std::set<int>> adj;

public:

    Module() = delete;

    const std::map<int, std::set<int>> &getAdj() const;

    Module(const std::string &name, Level level, int maxNumVertices);

    base::dynamic_bitset<> accessioned_entity_vertices;

    void addVertex(int vertex);

    void addVertex(int interactomeIndex, unsigned int moduleBitsetIndex);

    void addEdge(int index1, int index2);

    void addEdges(std::vector<std::pair<int, int>> edges);

    bool hasVertex(int index);

    std::vector<int> getVertices();

    std::set<int> getNeighbors(int index);

    const std::string &getName() const;

    Level getLevel() const;
};


#endif //PROTEOFORMNETWORKS_MODULE_HPP
