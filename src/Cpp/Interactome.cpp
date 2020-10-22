#include "Interactome.hpp"

Level Interactome::get_type(int index) {
    if (index <= end_indexes[genes])
        return genes;
    else if (index <= end_indexes[proteins])
        return proteins;
    else if (index <= end_indexes[proteoforms])
        return proteoforms;
    else
        return SimpleEntity;
}

Level Interactome::get_type(std::string_view name) {
    return get_type(vertices.index(name.data()));
}

bool Interactome::isGene(int index) {
    return index <= end_indexes[genes];
}

bool Interactome::isProtein(int index) {
    return start_indexes[proteins] <= index && index <= end_indexes[proteins];
}

bool Interactome::isProteoform(int index) {
    return start_indexes[proteoforms] <= index && index <= end_indexes[proteoforms];
}

bool Interactome::isSimpleEntity(int index) {
    return index >= start_indexes[SimpleEntity];
}

void Interactome::readEdges(std::string_view file_edges) {
    std::cout << "Reading edges...\n";
    std::ifstream f;
    f.open(file_edges.data());

    int i1, i2;
    while (f >> i1 >> i2) {
        if (adj.find(i1) == adj.end()) {
            std::set<int> s = {i2};
            adj.emplace(i1, s);
        } else {
            adj[i1].insert(i2);
        }

        if (adj.find(i2) == adj.end()) {
            std::set<int> s = {i1};
            adj.emplace(i2, s);
        } else {
            adj[i2].insert(i1);
        }
    }
    std::cout << "Complete.\n\n";
}

void Interactome::readRanges(std::string_view file_ranges) {
    std::cout << "Reading index ranges.\n";
    std::ifstream f;
    f.open(file_ranges.data());

    int start_index, end_index;
    while (f >> start_index >> end_index) {
        start_indexes.push_back(start_index);
        end_indexes.push_back(end_index);
    }
    std::cout << "Complete.\n\n";
}

int Interactome::index(std::string_view name) {
    return vertices.index(name.data());
}


void Interactome::readGenesToProteins(std::string_view file_proteins_to_genes) {
    std::cout << "Reading genes to proteins...\n";
    std::ifstream f(file_proteins_to_genes.data());

    if (!f.is_open()) {
        std::string message = "Error reading file ";
        message += file_proteins_to_genes;
        message += " at: ";
        message += __FUNCTION__;
        throw std::runtime_error(message);
    }

    std::string protein_name, gene_name;
    while (f >> protein_name >> gene_name) {
        int gene_index = index(gene_name), protein_index = index(protein_name);
        if (genes_to_proteins.find(gene_index) == genes_to_proteins.end()) {
            genes_to_proteins.emplace(gene_index, std::vector<int>(1, protein_index));
        } else {
            genes_to_proteins[gene_index].push_back(protein_index);
        }
    }
    std::cout << "Complete.\n\n";
}

void Interactome::readProteinsToProteoforms(std::string_view file_proteins_to_proteoforms) {
    std::cout << "Reading proteins to proteoforms.\n";
    std::ifstream f(file_proteins_to_proteoforms.data());

    if (!f.is_open()) {
        std::string message = "Error reading file ";
        message += file_proteins_to_proteoforms;
        message += " at: ";
        message += __FUNCTION__;
        throw std::runtime_error(message);
    }

    std::string protein_name, proteoform_name;
    while (f >> protein_name >> proteoform_name) {
        int proteoform_index = index(proteoform_name), protein_index = index(protein_name);
        if (proteins_to_proteoform.find(protein_index) == proteins_to_proteoform.end()) {
            proteins_to_proteoform.emplace(protein_index, proteoform_index);
        } else {
            proteins_to_proteoform[protein_index].push_back(proteoform_index);
        }
    }
    std::cout << "Complete.\n\n";
}

Interactome::Interactome(std::string_view file_vertices, std::string_view file_edges, std::string_view file_ranges,
                         std::string_view file_proteins_to_genes, std::string_view file_proteins_to_proteoforms)
        : vertices(file_vertices) {
    readRanges(file_ranges);
    readEdges(file_edges);
    readGenesToProteins(file_proteins_to_genes);
    readProteinsToProteoforms(file_proteins_to_proteoforms);
    std::cout << "Finished reading interactome.\n\n";
}

std::vector<int> Interactome::getProteins(int gene_index) {
    return genes_to_proteins[gene_index];
}

std::vector<int> Interactome::getProteoforms(int protein_index) {
    return proteins_to_proteoform[protein_index];
}

// Returns all types of neighbors of the vertex
std::set<int> Interactome::getNeighbors(int index) {
    return adj[index];
}

std::set<int> Interactome::getProteinNeighbors(int index) {
    std::set<int> neighbors;
    for (auto neighbor : adj[index]) {
        if (isProtein(neighbor)) neighbors.insert(neighbor);
    }
    return neighbors;
}

std::set<int> Interactome::getProteoformNeighbors(int index) {
    std::set<int> neighbors;
    for (auto neighbor : adj[index]) {
        if (isProteoform(neighbor)) neighbors.insert(neighbor);
    }
    return neighbors;
}

std::set<int> Interactome::getSimpleEntityNeighbors(int index) {
    std::set<int> neighbors;
    for (auto neighbor : adj[index]) {
        if (isSimpleEntity(neighbor)) neighbors.insert(neighbor);
    }
    return neighbors;
}

bool has(const std::vector<int> &indexes, int index) {
    auto lower = lower_bound(indexes.begin(), indexes.end(), index);
    return lower != indexes.end() && *lower == index;
}

// Returns edges between the selected vertices
std::vector<std::pair<int, int>> Interactome::getInteractions(std::vector<int> indexes) {
    std::vector<std::pair<int, int>> interactions;

    std::sort(indexes.begin(), indexes.end());

    for (int index : indexes) {
        for (int neighbor : getNeighbors(index)) {
            if (has(indexes, neighbor)) {
                interactions.push_back(std::make_pair(index, neighbor));
            }
        }
    }

    return interactions;
}

int Interactome::getNumVertices(){
    return vertices.size();
}

int Interactome::getStartIndexGenes() {
    return start_indexes[genes];
}

int Interactome::getEndIndexGenes() {
    return end_indexes[genes];
}

int Interactome::getStartIndexProteins() {
    return start_indexes[proteins];
}

int Interactome::getEndIndexProteins(){
    return end_indexes[proteins];
}

int Interactome::getStartIndexProteoforms(){
    return start_indexes[proteoforms];
}

int Interactome::getEndIndexProteoforms(){
    return end_indexes[proteoforms];
}