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

	vs getModifications(std::string proteoform) {
		vs modifications;
		std::sregex_iterator it(proteoform.begin(), proteoform.end(), RGX_MODIFICATION);
		std::sregex_iterator end;
		while (it != end) {
			if (it->str().find(';') || it->str().find(',')) {
				modifications.push_back(it->str().substr(1));
			}
			else {
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

	vs readProteoforms(const std::string& path_file_mapping, bool hasHeader) {
		std::ifstream map_file(path_file_mapping.data());
		std::string entity, leftover;
		uss temp_set;
		vs index_to_entities;

		if (!map_file.is_open()) {
			throw std::runtime_error("Could not open file " + path_file_mapping);
		}

		if (hasHeader) {
			getline(map_file, leftover);  // Skip header line
		}
		while (map_file.peek() != EOF) {
			if (map_file.peek() == '[' || map_file.peek() == '\"') {
				entity = proteoform::readProteoformFromNeo4jCsv(map_file);
			}
			else {
				getline(map_file, entity, ',');  // Read entity
			}
			getline(map_file, leftover);  // Read rest of line
			temp_set.insert(entity);
		}
		index_to_entities = convert_uss_to_vs(temp_set);

		return index_to_entities;
	}

	const measures_result calculateModificationsPerProteoform(const vs& proteoforms) {
		double min = 1.0;
		double max = 0.0;
		double avg = 0.0;
		double sum = 0.0;
		for (const auto& proteoform : proteoforms) {
			double num = static_cast<double>(proteoform::getModifications(proteoform).size());
			if (num < min)
				min = num;
			if (num > max)
				max = num;
			sum += num;
		}
		avg = sum / proteoforms.size();
		return { min, max, avg };
	}

}  // namespace proteoform