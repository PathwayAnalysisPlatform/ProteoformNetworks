#ifndef BIMAP_HPP_
#define BIMAP_HPP_

#include <fstream>
#include <bitset>
#include <string_view>
#include <cstring>
#include <set>
#include "types.hpp"

std::string rtrim(std::string &s);

vs convert_uss_to_vs(const uss &a_set);

vs createIntToStr(std::string_view path_file, bool has_header = true, int column_index = 0, int total_num_columns = 1);

umsi createStrToInt(const vs &index_to_entities);

class Bimap_str_int {

    umsi stoi; // string to int

    vs itos; // int to string

public:

    // Create from list, without header, list comes already sorted, there is only one column in the file.
    Bimap_str_int(std::string_view file_elements);

    // Create bimap from element to index in a unique sorted array.
    // Creates a bimap of the elements in the selected column.
    // Column index starts counting at 0
    Bimap_str_int(std::string_view file_elements, bool has_header, int column_index, int total_num_columns);

    Bimap_str_int(const vs &index_to_entities);

    int size() const { return itos.size(); }

    bool has(const std::string &key) const { return stoi.find(key) != stoi.end(); };

    int index(const std::string &key) const { return stoi.at(key); };

    std::string name(const int index) const { return itos[index]; };

};

#endif // !BIMAP_HPP_
