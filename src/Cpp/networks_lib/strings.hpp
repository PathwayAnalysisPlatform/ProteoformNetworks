#ifndef PROTEOFORMNETWORKS_STRINGS_HPP
#define PROTEOFORMNETWORKS_STRINGS_HPP

#include <string>
#include <cctype>

std::string rtrim(std::string &s){
    if(s.length() > 0){
        auto it = s.end()-1;
        while(isspace(*it)){
            s.erase(it);
            it = s.end()-1;
        }
    }
    return s;
}

#endif //PROTEOFORMNETWORKS_STRINGS_HPP
