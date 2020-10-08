#include <iostream>
#include "networks.hpp"

// Calculates and create a file with the sizes of all trait All_modules at all levels (genes, proteins or proteoforms)
std::map<const std::string, um<int, int>> calculate_and_report_sizes(std::string_view path_reports,
                                                                     const std::map<const char *, const All_modules> &all_modules,
                                                                     const bimap_str_int &groups) {
    std::map<const std::string, um<int, int>> sizes;

    for (const auto &level : config::LEVELS) {
        std::string file_name = path_reports.data() + static_cast<std::string>("module_sizes_") + level + ".tsv";
        std::ofstream output(file_name);

        if (!output.is_open()) {
            std::string message = "Cannot open report file at ";
            std::string function = __FUNCTION__;
            throw std::runtime_error(message + function);
        }

        std::cerr << "Writing file: " << file_name << std::endl;

        output << "MODULE\tSIZE\n";
        for (int group_index = 0; group_index < all_modules.at(level).group_to_members.size(); group_index++) {
            output << groups.int_to_str[group_index] << "\t"
                   << all_modules.at(level).group_to_members[group_index].count()
                   << "\n";
            sizes[level][group_index] = all_modules.at(level).group_to_members[group_index].count();
        }
        output.close();
    }

    return sizes;
}
