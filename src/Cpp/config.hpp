#ifndef PROTEOFORMNETWORKS_CONFIG_HPP
#define PROTEOFORMNETWORKS_CONFIG_HPP

#include <string>

namespace config {
    const std::string MODULES_FILE_END = "_modules.tsv";

    enum levels {
        genes, proteins, proteoforms
    };
    const std::vector<const char *> LEVELS{"genes",
                                           "proteins",
                                           "proteoforms"};
}


#endif //PROTEOFORMNETWORKS_CONFIG_HPP
