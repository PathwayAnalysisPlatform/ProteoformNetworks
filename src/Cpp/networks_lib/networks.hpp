#ifndef PROTEOFORMNETWORKS_NETWORKS_HPP
#define PROTEOFORMNETWORKS_NETWORKS_HPP

#include <string_view>
#include "types.hpp"
#include "bimap_str_int.hpp"

ummii loadInteractionNetwork(std::string_view path_file_interactions,
                             const bimap_str_int &entities,
                             bool has_header_row);

#endif //PROTEOFORMNETWORKS_NETWORKS_HPP
