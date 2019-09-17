#include "gtest/gtest.h"
#include "phegeni.hpp"

class PhegeniCreatePheGenIProteinModulesFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        // Create gene modules with the slice of PheGenI and the slice of Reactome genes
//        trait_modules gene_modules = loadPheGenIGeneModules(
//                path_file_phegeni, path_file_genes);
//        trait_modules protein_modules = createPheGenIProteinModules(gene_modules, path_file_genes, );
    }
    std::string path_file_phegeni = "../../../../resources/PheGenI/PheGenI_Association_genome_wide_significant_slice.txt";
    std::string path_file_genes = "../../../../resources/Reactome/v70/genes_slice.csv";
};

// Check the number of modules is correct
TEST(PhegeniCreatePheGenIProteinModulesFixture, CorrectModules) {

}
