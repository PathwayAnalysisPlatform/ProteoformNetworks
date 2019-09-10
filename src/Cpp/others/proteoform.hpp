#ifndef PROTEOFORM_H
#define PROTEOFORM_H

#include <string>
#include <regex>
#include <set>
#include <bitset>
#include <fstream>
#include <scores.hpp>

#include "bimap_str_int.hpp"
#include "types.hpp"
#include "conversions.hpp"

namespace proteoform {

// Method to get convert from Neo4j csv format to our simple custom format

// Method to get the accesion from a proteoform string in simple format

const std::regex RGX_ACCESSION_DELIMITER{"[;-]"};
const std::regex RGX_MODIFICATION{"[;,]\\d{5}"};

std::string getAccession(const std::string& proteoform);

// Method to verify if a proteoform has modifications
bool isModified(const std::string& proteoform);

template <size_t S>
std::bitset<S> getSetOfModifiedProteoforms(const bimap_str_int& proteoforms) {
	std::bitset<S> modified_proteoforms;

	for (int I = 0; I < proteoforms.int_to_str.size(); I++) {
		if (isModified(proteoforms.int_to_str[I])) {
			modified_proteoforms.set(I);
		}
	}

	return modified_proteoforms;
}

vs getModifications(std::string proteoform);

std::string readProteoformFromNeo4jCsv(std::ifstream& fs);

const measures_result calculateModificationsPerProteoform(const vs& proteoforms);

} // namespace proteoform

#endif /* PROTEOFORM_H */

