#ifndef PROTEOFORMNETWORKS_NETWORKS_HPP
#define PROTEOFORMNETWORKS_NETWORKS_HPP

#include <string_view>
#include "types.hpp"
#include "bimap_str_int.hpp"

vusi loadInteractionNetwork(std::string_view path_file_interactions,
                            const bimap_str_int &entities,
                            bool has_header_row = true);

struct load_modules_result {
    modules entity_modules;
    bimap_str_int groups;
    bimap_str_int members;
};

load_modules_result loadModules(std::string_view path_file_modules, bool has_header = true);

void writeModules(std::string_view path_file_modules,
                  const modules &entity_modules,
                  const bimap_str_int &groups,
                  const bimap_str_int &members);

modules removeDisconnectedMembers(modules &modules, const bimap_str_int &groups, const bimap_str_int &members,
                                  const vusi &interactions);

#endif //PROTEOFORMNETWORKS_NETWORKS_HPP
