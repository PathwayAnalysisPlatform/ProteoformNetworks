#include "Module.hpp"

Module::Module() {
    throw std::runtime_error("This constructor should not be called.");
}

Module::Module(const std::string &name, Level level) : name(name), level(level) {}

void Module::addVertex(int vertex) {
    adj[vertex];
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



