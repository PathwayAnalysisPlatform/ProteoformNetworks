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

std::vector<std::string> getModifications(std::string proteoform) {
   std::vector<std::string> modifications;
   std::sregex_iterator it(proteoform.begin(), proteoform.end(), RGX_MODIFICATION);
   std::sregex_iterator end;
   while (it != end) {
      if (it->str().find(';') || it->str().find(',')) {
         modifications.push_back(it->str().substr(1));
      } else {
         modifications.push_back((*it)[0]);
      }
      it++;
   }
   return modifications;
}

std::string readProteoformFromNeo4jCsv(std::ifstream& fs) {
   std::string proteoform;
   if (fs.peek() == '\"')  // Read initial "
      fs.get();
   fs.get();  // Read initial " or [
   getline(fs, proteoform, ']');

   // Remove the extra ',' in the proteoform representation comming directly from Neo4j queries
   std::size_t index = proteoform.find_first_of(',');
   if (index != std::string::npos)
      proteoform.erase(index, 1);
   return proteoform;
}

}  // namespace proteoform