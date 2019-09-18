#include "gtest/gtest.h"
#include "phegeni.hpp"

class PhegeniCreatePheGenIProteinModulesFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        // Create gene modules with the slice of PheGenI and the slice of Reactome genes
        genes = createBimap(path_file_genes);
        proteins = createBimap(path_file_proteins);
        auto ret = loadPheGenIGenesAndTraits(path_file_phegeni, genes);
        traits = ret.groups;
        modules gene_modules = loadPheGenIGeneModules(path_file_phegeni, genes, traits);
        modules protein_modules = createPheGenIProteinModules(gene_modules,
                                                              genes,
                                                              proteins,
                                                              traits,
                                                              path_file_mapping,
                                                              path_file_protein_interactions);
    }

    std::string path_file_phegeni = "../../../../resources/PheGenI/PheGenI_Association_genome_wide_significant_slice.txt";
    std::string path_file_genes = "../../../../resources/Reactome/v70/Genes/genes_slice.csv";
    std::string path_file_proteins = "../../../../resources/Reactome/v70/Proteins/proteins_slice.csv";
    std::string path_file_protein_interactions = "../../../../Reactome/v70/Proteins/proteinInternalEdges_slice.csv";
    std::string path_file_mapping = "../../../../resources/UniProt/mapping_protein_to_genes_slice.tab";
    bimap_str_int genes, proteins, traits;
};

// Check the number of modules is correct
TEST(PhegeniCreatePheGenIProteinModulesFixture, CorrectModules) {
    //
}
