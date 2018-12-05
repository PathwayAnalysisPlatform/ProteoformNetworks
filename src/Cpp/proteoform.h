#ifndef PROTEOFORM_H
#define PROTEOFORM_H

#include <string>
#include <regex>
#include <set>
#include <bitset>

namespace proteoform {

// Method to get convert from Neo4j csv format to our simple custom format

// Method to get the accesion from a proteoform string in simple format

const std::regex RGX_ACCESSION_DELIMITER{"[;-]"};
const std::regex RGX_MODIFICATION{"[;,]\\d{5}"};

std::string getAccession(const std::string& proteoform) {
   std::smatch match_end_of_accession;
   if (!regex_search(proteoform, match_end_of_accession, RGX_ACCESSION_DELIMITER)) {
      return proteoform;
   }
   return proteoform.substr(0, match_end_of_accession.position(0));
}

// Method to verify if a proteoform has modifications
bool isModified(const std::string& proteoform) {
   std::smatch modification;
   return regex_search(proteoform, modification, RGX_MODIFICATION);
}

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

std::set<std::string> getModifications(std::string proteoform) {
   std::set<std::string> modifications;
   std::sregex_iterator it(proteoform.begin(), proteoform.end(), RGX_MODIFICATION);
   std::sregex_iterator end;
   while (it != end) {
      if (it->str().find(';') || it->str().find(',')) {
         modifications.insert(it->str().substr(1));
      } else {
         modifications.insert((*it)[0]);
      }
      it++;
   }
   return modifications;
}

} // namespace proteoform

#endif /* PROTEOFORM_H */

