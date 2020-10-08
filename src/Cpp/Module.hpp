#ifndef PROTEOFORMNETWORKS_MODULE_HPP
#define PROTEOFORMNETWORKS_MODULE_HPP


#include <string>
#include <map>
#include <set>
#include <vector>
#include "Interactome.hpp"

class Module {
    std::string name;
    Level level;
    std::map<int, std::set<int>> adj;


public:

    Module();

    Module(const std::string &name, Level level);

    void addVertex(int vertex);

    void addEdge(int index1, int index2);

    void addEdges(std::vector<std::pair<int, int>> edges);

    bool hasVertex(int index);

    std::vector<int> getVertices();

    std::set<int> getNeighbors(int index);

    const std::string &getName() const;

    Level getLevel() const;
};


#endif //PROTEOFORMNETWORKS_MODULE_HPP
