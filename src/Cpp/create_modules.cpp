#include <iostream>
#include "overlap_analysis.hpp"
#include "Interactome.hpp"

int main(int argc, char *argv[]) try {

    if (argc < 4) {
        std::cerr << "Missing arguments. Expected: 4 arguments:\n\n"
                  << " * - [1] File with Gene sets to create disease modules\n"
                  << " * - [2] File of Interactome edges with indexed vertices\n"
                  << " * - [3] File with start and end vertex indexes\n"
                  << " * - [4] Mapping from proteins to genes\n"
                  << " * - [5] Mapping from proteins to proteoforms\n"
                  << " * - [6] Output path";
        throw std::runtime_error("Missing arguments.");
        return 0;
    }

    std::string file_phegeni = argv[1];
    std::string file_interactome = argv[2];
    std::string file_indexes = argv[3];
    std::string file_map_proteins_to_genes = argv[4];
    std::string file_map_proteins_to_proteoforms = argv[5];
    std::string output_path = argv[6];

    Interactome interactome = new Interactome(file_interactome, file_indexes);

    // Create or read module files at the three levels: all in one, and single module files.
    std::map<const char *, const All_modules> all_modules = get_or_create_modules(
            argv[1],
            argv[2],
            argv[3],
            argv[4],
            argv[5],
            argv[6],
            argv[7],
            argv[8],
            argv[9],
            argv[10],
            argv[11],
            true);
}
catch (const std::exception &ex) {
    std::cout << ex.what() << "\n";
}