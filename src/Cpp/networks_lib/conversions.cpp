#include "conversions.hpp"

// The result vector is sorted lexicographycally
vs convert_uss_to_vs(const uss& a_set) {
    vs result;
    result.assign(a_set.begin(), a_set.end());
    return result;
}

entity_mapping readMapping(std::string_view path_file_mapping) {
    ummss source_to_destinations;
    ummss destination_to_sources;
//    ifstream file_mapping(path_file_mapping.data());
//    string protein, proteoform, leftover;
//
//    if (!file_mapping.is_open()) {
//        std::string message = "Error reading file at: ";
//        std::string function = __FUNCTION__;
//        throw runtime_error(message + function);
//    }
//
//    getline(file_mapping, protein);  // Discard the header line.
//    while (getline(file_mapping, proteoform, '\t')) {
//        getline(file_mapping, protein, '\t');
//        getline(file_mapping, leftover);
//        protein_to_proteforms.emplace(protein, proteoform);
//        proteoform_to_proteins.emplace(proteoform, protein);
//    }
//
//    cerr << "Mapping proteins proteoforms loaded: " << protein_to_proteforms.size() << " = " << proteoform_to_proteins.size() << "\n";
    return { source_to_destinations, destination_to_sources };
}

