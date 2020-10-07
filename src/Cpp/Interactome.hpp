#ifndef PROTEOFORMNETWORKS_INTERACTOME_HPP
#define PROTEOFORMNETWORKS_INTERACTOME_HPP

#include <string_view>
#include <map>
#include <vector>
#include <set>
#include "config.hpp"

enum level {
    genes, proteins, proteoforms, SimpleEntity
};
std::vector<std::string> levels_str = {"genes", "proteins", "proteoforms", "SimpleEntity"};

class Interactome {
    std::map<int, std::set<int>> adj;

    std::vector<int> start_indexes;
    std::vector<int> end_indexes;

    int index(std::string_view name) {

    }

    int name(int index) {

    }

    level get_type(int index);
    level get_type(std::string_view name);

    void read_interactions(std::string_view file_interactions);
    void read_indexes(std::string_view file_indexes);

public:
    Interactome(std::string_view file_interactions, std::string_view file_indexes);
};


#endif //PROTEOFORMNETWORKS_INTERACTOME_HPP
