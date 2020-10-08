#include <iostream>
#include "overlap_analysis.hpp"
#include "Interactome.hpp"
#include "Module.hpp"


void saveModules(std::map<std::string, Module> modules, const std::string &output_path);

// Create or read module files at the three levels: all in one, and single module files.
// One file with the list of diseases and for each disease a module file for each level
std::map<std::string, Module> createGeneModules(std::string_view file_phegeni,
                                                Interactome interactome,
                                                const std::string &path_output) {

    // Read traits and genes in Phegeni file. The genes used are the ones from Reactome, if the gene read is not there it is ignored.
    std::ifstream file_phegen(file_phegeni.data());
    std::string line, field, trait, gene_name, gene_name_2;
    std::string p_value_str;
    long double p_value;

    std::map<std::string, Module> modules;

    if (!file_phegen.is_open()) {
        std::string message = "Cannot open path_file_phegeni at ";
        std::string function = __FUNCTION__;
        throw std::runtime_error(message + function);
    }

    getline(file_phegen, line);                  // Read header line
    while (getline(file_phegen, field, '\t')) {  // Read #
        getline(file_phegen, trait, '\t');        // Read Trait
        getline(file_phegen, field, '\t');        // Read SNP rs
        getline(file_phegen, field, '\t');        // Read Context
        getline(file_phegen, gene_name, '\t');         //	Gene
        getline(file_phegen, field, '\t');        //	Gene ID
        getline(file_phegen, gene_name_2, '\t');        //	Gene 2
        getline(file_phegen, field, '\t');        //	Gene ID 2
        getline(file_phegen, field, '\t');        // Read Chromosome
        getline(file_phegen, field, '\t');        // Read Location
        getline(file_phegen, p_value_str, '\t');  // Read P-Value
        getline(file_phegen,
                line);               // Skip header line leftoever: Source,	PubMed,	Analysis ID,	Study ID,	Study Name

        if (modules.find(trait) == modules.end()) {
            modules.emplace(trait, Module(trait, genes));
        }
        modules[trait].addVertex(interactome.index(gene_name));
        for (int simpleEntity : interactome.getSimpleEntityNeighbors(interactome.index(gene_name)))
            modules[trait].addVertex(simpleEntity);

        modules[trait].addVertex(interactome.index(gene_name_2));
        for (int simpleEntity : interactome.getSimpleEntityNeighbors(interactome.index(gene_name_2)))
            modules[trait].addVertex(simpleEntity);
    }

    modules[trait].addEdges(interactome.getInteractions(modules[trait].getVertices()));

    for (auto entry : modules) {
        Module &module = entry.second;
        std::ofstream f(path_output.data() + module.getName() + "_" + "genes.tsv");

        for (auto vertex : module.getVertices()) {
            if (module.getNeighbors(vertex).size() == 0) {
                f << vertex << "\t" << vertex << "\n";
            } else {
                for (auto neighbor : module.getNeighbors(vertex)) {
                    f << vertex << "\t" << neighbor << "\n";
                }
            }
        }
    }

    saveModules(modules, path_output);

    std::cerr << "Created " << modules.size() << " disease modules\n";

    return modules;
}

std::map<std::string, Module> createProteinModules(std::map<std::string, Module> gene_modules,
                                                   Interactome interactome,
                                                   const std::string &output_path) {

    std::map<std::string, Module> protein_modules;
    for (auto &entry : gene_modules) {
        Module &gene_module = entry.second;
        Module protein_module = Module(gene_module.getName(), proteins);

        std::set<int> candidate_proteins;
        for (auto vertex : gene_module.getVertices()) { // Gather all candidate proteins by converting ids
            if (interactome.isGene(vertex)) {
                for (auto protein : interactome.getProteins(vertex))
                    candidate_proteins.insert(protein);
            }
        }

        for (auto candidate_protein : candidate_proteins) { // Keep only those with interactors in the same set
            for (auto neighbor : interactome.getProteinNeighbors(candidate_protein)) {
                if (candidate_proteins.find(neighbor) != candidate_proteins.end()) {
                    protein_module.addVertex(candidate_protein);
                    break;
                }
            }
        }

        protein_module.addEdges(interactome.getInteractions(protein_module.getVertices()));

        protein_modules.emplace(protein_module.getName(), protein_module);
    }

    saveModules(protein_modules, output_path);

    return protein_modules;
}

std::map<std::string, Module> createProteoformModules(std::map<std::string, Module> protein_modules,
                                                      Interactome interactome,
                                                      const std::string &output_path) {

    std::map<std::string, Module> proteoform_modules;
    for (auto &entry : protein_modules) {
        Module &protein_module = entry.second;
        Module proteoform_module = Module(protein_module.getName(), proteoforms);

        std::set<int> candidate_proteoforms;
        for (auto vertex : protein_module.getVertices()) { // Gather all candidate proteins by converting ids
            if (interactome.isProtein(vertex)) {
                for (auto proteoform : interactome.getProteoforms(vertex))
                    candidate_proteoforms.insert(proteoform);
            }
        }

        for (auto candidate_proteoform : candidate_proteoforms) { // Keep only those with interactors in the same set
            for (auto neighbor : interactome.getProteoformNeighbors(candidate_proteoform)) {
                if (candidate_proteoforms.find(neighbor) != candidate_proteoforms.end()) {
                    proteoform_module.addVertex(candidate_proteoform);
                    break;
                }
            }
        }

        proteoform_module.addEdges(interactome.getInteractions(proteoform_module.getVertices()));

        proteoform_modules.emplace(proteoform_module.getName(), proteoform_module);
    }

    saveModules(proteoform_modules, output_path);

    return proteoform_modules;

}

void saveModules(std::map<std::string, Module> modules, const std::string &output_path) {
    for (auto entry : modules) {
        Module &module = entry.second;
        std::ofstream f(output_path.data() + module.getName() + "_" + LEVELS[module.getLevel()] + ".tsv");

        for (auto vertex : module.getVertices()) {
            if (module.getNeighbors(vertex).size() == 0) {
                f << vertex << "\t" << vertex << "\n";
            } else {
                for (auto neighbor : module.getNeighbors(vertex)) {
                    f << vertex << "\t" << neighbor << "\n";
                }
            }
        }
    }
}

void createModules(std::string_view file_phegeni,
                   Interactome interactome,
                   const std::string &output_path) {
    std::map<std::string, Module> gene_modules = createGeneModules(file_phegeni, interactome, output_path);
    std::map<std::string, Module> protein_modules = createProteinModules(gene_modules, interactome, output_path);
    std::map<std::string, Module> proteoform_modules = createProteoformModules(protein_modules, interactome,
                                                                               output_path);
}


int main(int argc, char *argv[]) try {

    if (argc < 4) {
        std::cerr << "Missing arguments. Expected: 4 arguments:\n\n"
                  << " * - [1] File with Gene sets to create disease modules\n"
                  << " * - [2] File of Interactome vertices\n"
                  << " * - [3] File of Interactome edges (with indexed vertices)\n"
                  << " * - [4] File with start and end vertex indexes\n"
                  << " * - [5] Mapping from proteins to genes\n"
                  << " * - [6] Mapping from proteins to proteoforms\n"
                  << " * - [7] Output path";
        throw std::runtime_error("Missing arguments.");
        return 0;
    }

    std::string file_phegeni = argv[1];
    std::string file_vertices = argv[2];
    std::string file_edges = argv[3];
    std::string file_indexes = argv[4];
    std::string file_proteins_to_genes = argv[5];
    std::string file_proteins_to_proteoforms = argv[6];
    std::string output_path = argv[7];

    Interactome interactome(file_vertices, file_edges, file_indexes, file_proteins_to_genes,
                            file_proteins_to_proteoforms);

    createModules(file_phegeni, interactome, output_path);
}
catch (const std::exception &ex) {
    std::cout << ex.what() << "\n";
}