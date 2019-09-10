#ifndef ENTITY_HPP_
#define ENTITY_HPP_

#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string_view>

#include "proteoform.hpp"
#include "conversions.hpp"

enum entities {
	GENES,
	PROTEINS,
	PROTEOFORMS
};

#endif /* ENTITY_HPP_ */