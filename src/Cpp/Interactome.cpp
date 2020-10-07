#include <fstream>
#include "Interactome.hpp"

level Interactome::get_type(int index) {
    if (index <= end_indexes[genes])
        return genes;
    else if (index <= end_indexes[proteins])
        return proteins;
    else if (index <= end_indexes[proteoforms])
        return proteoforms;
    else
        return SimpleEntity;
}

level Interactome::get_type(std::string_view name) {
    return get_type(index(name));
}

void Interactome::read_indexes(std::string_view file_indexes) {
    std::ifstream f;
    f.open(file_indexes.data());

    int i1, i2;
    while(f >> i1 >> i2){
        if(adj.find(i1) == adj.end()){
            std::set<int> s = {i2};
            adj.emplace(i1, s);
        } else{
            adj[i1].insert(i2);
        }

        if(adj.find(i2) == adj.end()){
            std::set<int> s = {i1};
            adj.emplace(i2, s);
        } else{
            adj[i2].insert(i1);
        }
    }
}

void Interactome::read_interactions(std::string_view file_interactions){
    std::ifstream f;
    f.open(file_interactions.data());

    int start_index, end_index;
    while(f >> start_index >> end_index){
        start_indexes.push_back(start_index);
        end_indexes.push_back(end_index);
    }
}

Interactome::Interactome(std::string_view file_interactions, std::string_view file_indexes) {
    read_indexes(file_indexes);
    read_interactions(file_interactions);
}
