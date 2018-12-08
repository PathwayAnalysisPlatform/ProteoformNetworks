#include "proteoform.hpp"

namespace proteoform {

std::string getAccession(const std::string& proteoform) {
   std::smatch match_end_of_accession;
   if (!regex_search(proteoform, match_end_of_accession, RGX_ACCESSION_DELIMITER)) {
      return proteoform;
   }
   return proteoform.substr(0, match_end_of_accession.position(0));
}

bool isModified(const std::string& proteoform) {
   std::smatch modification;
   return regex_search(proteoform, modification, RGX_MODIFICATION);
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

}