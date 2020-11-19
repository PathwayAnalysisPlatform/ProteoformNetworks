#include "Module.hpp"

Module::Module(const std::string &name, Level level, int maxNumVertices) :
        name(name),
        level(level),
        accessioned_entity_vertices(base::dynamic_bitset<>(maxNumVertices)) {

    //    std::cout << "Module for " << LEVELS[level] << " constructed with " << maxNumVertices << std::endl;
}

// Add vertex to the adjacency list
void Module::addVertex(int vertex) {
    adj[vertex];
}

// Add vertex to the adjacency list and to the bitset for overlap operations
void Module::addVertex(int interactomeIndex, unsigned int moduleBitsetIndex) {
    adj[interactomeIndex];
    if(moduleBitsetIndex >= accessioned_entity_vertices.size()){
        std::cerr << "Tried to add entity out of index range:" << moduleBitsetIndex << ". Max index: " << accessioned_entity_vertices.size()-1 << std::endl;
        std::cerr << "Level: " << LEVELS[level] << std::endl;
    } else{
        accessioned_entity_vertices[moduleBitsetIndex] = true;
    }
}

void Module::addEdge(int index1, int index2) {
    addVertex(index1);
    addVertex(index2);

    if (index1 != index2) {
        adj[index1].insert(index2);
        adj[index2].insert(index1);
    }
}

void Module::addEdges(std::vector<std::pair<int, int>> edges) {
    for (auto edge : edges) {
        addEdge(edge.first, edge.second);
    }
}

bool Module::hasVertex(int index) {
    return adj.find(index) != adj.end();
}

std::vector<int> Module::getVertices() {
    std::vector<int> vertices;
    for (auto &entry :  adj)
        vertices.push_back(entry.first);
    return vertices;
}

std::set<int> Module::getNeighbors(int index) {
    return adj[index];
}

const std::string &Module::getName() const {
    return name;
}

Level Module::getLevel() const {
    return this->level;
}

const std::map<int, std::set<int>> &Module::getAdj() const {
    return adj;
}



