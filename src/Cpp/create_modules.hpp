#ifndef PROTEOFORMNETWORKS_CREATE_MODULES_HPP
#define PROTEOFORMNETWORKS_CREATE_MODULES_HPP

#include <string_view>
#include "Interactome.hpp"
#include "Module.hpp"


// Create or read module files at the three levels: all in one, and single module files.
std::map<std::string, Module> createGeneModules(std::string_view file_phegeni,
                                                Interactome interactome,
                                                const std::string &path_output);

std::map<std::string, Module> createProteinModules(std::map<std::string, Module> gene_modules,
                                                   Interactome interactome,
                                                   const std::string &path_output);

std::map<std::string, Module> createProteoformModules(std::map<std::string, Module> protein_modules,
                                                      Interactome interactome,
                                                      const std::string &output_path);

void createModules(std::string_view file_phegeni, Interactome interactome, const std::string &output_path);

void saveModules(std::map<std::string, Module> modules, const std::string &output_path);

#endif //PROTEOFORMNETWORKS_CREATE_MODULES_HPP
