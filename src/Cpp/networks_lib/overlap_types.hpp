#ifndef PROTEOFORMNETWORKS_OVERLAP_TYPES_HPP
#define PROTEOFORMNETWORKS_OVERLAP_TYPES_HPP

#include <types.hpp>
#include <bimap_str_int.hpp>

struct bidirectional_mapping {
    ummss first_to_second;
    ummss second_to_first;
};

struct modules {
    msb group_to_members;
    msb member_to_groups;
};

#endif //PROTEOFORMNETWORKS_OVERLAP_TYPES_HPP
