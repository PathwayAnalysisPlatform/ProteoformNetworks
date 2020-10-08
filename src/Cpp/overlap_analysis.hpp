#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <iostream>
#include <ctime>
#include <sys/stat.h>
#include <utility>

#include "scores.hpp"
#include "bimap_str_int.hpp"
#include "others/uniprot.hpp"

#include <windows.h>

inline bool file_exists(const std::string &name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}


#endif /* OVERLAP_H_ */