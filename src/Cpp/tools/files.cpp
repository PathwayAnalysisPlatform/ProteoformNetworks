#include <sys/stat.h>
#include "files.hpp"

namespace files{
    bool file_exists(const std::string &name) {
        struct stat buffer{};
        return (stat(name.c_str(), &buffer) == 0);
    }

    std::string GetExeFileName() {
        char buffer[MAX_PATH];
        GetModuleFileName(nullptr, buffer, MAX_PATH);
        return std::string(buffer);
    }

    std::string GetExePath() {
        std::string f = GetExeFileName();
        return f.substr(0, f.find_last_of("\\/"));
    }
}