#include <string>
#include <cstdio>
#include <iostream>
#include "overlap_analysis.hpp"


/*Performs the Overlap analysis calculations and creates files with the results.
 *
 * INPUT:
 * - [1] Modules file
 * - [2] Gene list (has header row)
 * - [3] Proteins list (has header row)
 * - [4] Proteoform list (does NOT have header row)
 * - [5] Gene interactions
 * - [6] Protein interactions
 * - [7] Proteoform interactions
 * - [8] Mapping from proteins to genes
 * - [9] Mapping from proteins to proteoforms
 * - [10] Output path
 *
 * */
int main(int argc, char *argv[]) try {
//	auto out = freopen("out.txt", "w", stdout);
//	auto err = freopen("err.txt", "w", stderr);

    if (argc < 11) {
        std::cerr << "Missing arguments. Expected: 10 arguments:\n\n"
                  << " * - [1] Modules file\n"
                  << " * - [2] Gene list\n"
                  << " * - [3] Proteins list\n"
                  << " * - [4] Proteoform list\n"
                  << " * - [5] Gene interactions\n"
                  << " * - [6] Protein interactions\n"
                  << " * - [7] Proteoform interactions\n"
                  << " * - [8] Mapping from proteins to genes\n"
                  << " * - [9] Mapping from proteins to proteoforms\n"
                  << " * - [10] Output path";
        throw std::runtime_error("Missing arguments.");
        return 0;
    }

    doOverlapAnalysis(argv[1],
                      argv[2],
                      argv[3],
                      argv[4],
                      argv[5],
                      argv[6],
                      argv[7],
                      argv[8],
                      argv[9],
                      argv[10]);

}
catch (const std::exception &ex) {
    std::cout << ex.what() << "\n";
}
