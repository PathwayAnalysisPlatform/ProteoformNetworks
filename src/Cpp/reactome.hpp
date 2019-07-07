#ifndef REACTOME_HPP_
#define REACTOME_HPP_

#include <unordered_map>
#include <bitset>
#include <cstdio>

#include "types.hpp"
#include "bimap.hpp"
#include "entity.hpp"

// Reactome v69
const size_t REACTOME_GENES = 24109;
const size_t REACTOME_PROTEINS = 10833;
const size_t REACTOME_PROTEOFORMS = 13991;

bimap loadReactomeEntities(std::string_view path_search_file);
vs loadReactomeEntitiesIndexToEntitites(std::string_view path_search_file);

#endif // !REACTOME_HPP_

