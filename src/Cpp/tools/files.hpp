//
// Created by luisp on 11/10/2020.
//

#ifndef PROTEOFORMNETWORKS_FILES_HPP
#define PROTEOFORMNETWORKS_FILES_HPP

#include <string>
#include <windows.h>

namespace files{
    bool file_exists(const std::string &name);

    std::string GetExeFileName();

    std::string GetExePath();
}

#endif //PROTEOFORMNETWORKS_FILES_HPP
