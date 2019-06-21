#ifndef PROTEOFORM_H
#define PROTEOFORM_H

#include <string>
#include <regex>
#include <set>
#include <bitset>
#include <fstream>

namespace proteoform {

// Method to get convert from Neo4j csv format to our simple custom format

// Method to get the accesion from a proteoform string in simple format

const std::regex RGX_ACCESSION_DELIMITER{"[;-]"};
const std::regex RGX_MODIFICATION{"[;,]\\d{5}"};

std::string getAccession(const std::string& proteoform);

// Method to verify if a proteoform has modifications
bool isModified(const std::string& proteoform);

template <size_t S>
std::bitset<S> getSetOfModifiedProteoforms(const std::vector<std::string>& proteoforms) {
   std::bitset<S> modified_proteoforms;

   for (int I = 0; I < proteoforms.size(); I++) {
      if (isModified(proteoforms[I])) {
         modified_proteoforms.set(I);
      }
   }

   return modified_proteoforms;
}

std::vector<std::string> getModifications(std::string proteoform);

std::string readProteoformFromNeo4jCsv(std::ifstream& fs);

} // namespace proteoform

#endif /* PROTEOFORM_H */

