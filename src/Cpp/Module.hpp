#ifndef PROTEOFORMNETWORKS_MODULE_HPP
#define PROTEOFORMNETWORKS_MODULE_HPP


#include <string>
#include <map>
#include <set>

class Module {
    std::string level;
    std::map<int, std::set<int>> adj;

    Module(std::string_view level) {

    }
};


#endif //PROTEOFORMNETWORKS_MODULE_HPP
