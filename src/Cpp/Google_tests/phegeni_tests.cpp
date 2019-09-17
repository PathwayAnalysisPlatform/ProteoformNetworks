#include "gtest/gtest.h"
#include <bimap_str_int.hpp>
#include "phegeni.hpp"

class PhegeniLoadPheGenISetsFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        genes = createBimap(path_file_genes);
        auto ret = loadPheGenIGenesAndTraits(path_file_phegeni, genes);
        phegeni_genes = ret.phegeni_genes;
        phegeni_traits = ret.phegeni_traits;
        modules = loadPheGenIGeneModules(path_file_phegeni, path_file_genes);
    }

    std::string path_file_phegeni = "../../../../resources/PheGenI/PheGenI_Association_genome_wide_significant_slice.txt";
    std::string path_file_genes = "../../../../resources/Reactome/v70/Genes/genes_slice.csv";
    bimap_str_int genes;
    bimap_str_int phegeni_genes;
    bimap_str_int phegeni_traits;
    trait_modules modules;
};

// Throws runtime error if file PheGenI is not found
TEST_F(PhegeniLoadPheGenISetsFixture, UnexistentFileRaisesError) {
    try {
        auto genesAndTraits = loadPheGenIGenesAndTraits("wrong_path", createBimap(path_file_genes));
        FAIL() << "Expected std::runtime_error";
    }
    catch (std::runtime_error const &err) {
        EXPECT_EQ(err.what(), std::string("Cannot open path_file_phegeni at loadPheGenIGenesAndTraits"));
    }
    catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}
// Check that genes not in Gene list are not added

TEST_F(PhegeniLoadPheGenISetsFixture, LoadPheGenIGenesAndTraits) {

    // Check that the traits are correct
    EXPECT_EQ(8, phegeni_traits.str_to_int.size());
    EXPECT_EQ(8, phegeni_traits.int_to_str.size());

    // Check that the genes are correct
    EXPECT_EQ(25, phegeni_genes.int_to_str.size());
    EXPECT_EQ(25, phegeni_genes.int_to_str.size());

    // Check genes are sorted
    EXPECT_EQ(1, phegeni_traits.str_to_int.at("\"Cholesterol, LDL\""));
    EXPECT_EQ("\"Cholesterol, LDL\"", phegeni_traits.int_to_str[1]);
    EXPECT_EQ(2, phegeni_traits.str_to_int.at("Bilirubin"));
    EXPECT_EQ("Bilirubin", phegeni_traits.int_to_str[2]);
    EXPECT_EQ(7, phegeni_traits.str_to_int.at("Wet Macular Degeneration"));
    EXPECT_EQ("Wet Macular Degeneration", phegeni_traits.int_to_str[7]);
}

// Test the trait to gene bitsets contains the right number of trait keys
// Test the gene to trait bitsets contains the right number of gene keys
TEST_F(PhegeniLoadPheGenISetsFixture, TraitKeys) {
    EXPECT_EQ(8, modules.traits_to_entities.size()) << "It has the wrong number of Traits as keys of the map.";
    EXPECT_EQ(25, modules.entities_to_traits.size()) << "It has the wrong number of GENES as keys of the map.";
}

// Test the trait to gene bitset contains the correct number of bits in the bitsets
// Test the gene to trait bitsets contain the correct number of bits in the bitsets
TEST_F(PhegeniLoadPheGenISetsFixture, BitsetSizes) {
    EXPECT_EQ(25, modules.traits_to_entities["Bilirubin"].size()) << "The bitsets size are not the number of genes.";
    EXPECT_EQ(8, modules.entities_to_traits["UGT1A1"].size()) << "The bitsets size are not the number of traits.";
}

// Check a trait has the right number of genes as members
TEST_F(PhegeniLoadPheGenISetsFixture, CorrectGeneMembers) {
    EXPECT_EQ(9, modules.traits_to_entities["Bilirubin"].count());
    EXPECT_TRUE(modules.traits_to_entities["Bilirubin"][genes.str_to_int["UGT1A1"]]); // From column GENE ID
    EXPECT_TRUE(modules.traits_to_entities["Bilirubin"][genes.str_to_int["UGT1A4"]]); // From column GENE ID 2
    EXPECT_FALSE(modules.traits_to_entities["Bilirubin"][genes.str_to_int["AAAA"]]);  // From column GENE ID
    EXPECT_FALSE(modules.traits_to_entities["Bilirubin"][genes.str_to_int["DDDD"]]);  // From column GENE ID 2

    EXPECT_EQ(2, modules.traits_to_entities["\"Cholesterol, HDL\""].count());
    EXPECT_TRUE(modules.traits_to_entities["\"Cholesterol, HDL\""][genes.str_to_int["HERPUD1"]]); // From column GENE ID
    EXPECT_TRUE(modules.traits_to_entities["\"Cholesterol, HDL\""][genes.str_to_int["CETP"]]); // From column GENE ID 2
    EXPECT_FALSE(modules.traits_to_entities["\"Cholesterol, HDL\""][genes.str_to_int["APOE"]]);

    EXPECT_EQ(9, modules.traits_to_entities["Bilirubin"].count());
    EXPECT_TRUE(modules.traits_to_entities["Bilirubin"][genes.str_to_int["UGT1A1"]]); // From column GENE ID
    EXPECT_TRUE(modules.traits_to_entities["Bilirubin"][genes.str_to_int["UGT1A4"]]); // From column GENE ID 2
    EXPECT_FALSE(modules.traits_to_entities["Bilirubin"][genes.str_to_int["AAAA"]]);  // From column GENE ID
    EXPECT_FALSE(modules.traits_to_entities["Bilirubin"][genes.str_to_int["DDDD"]]);  // From column GENE ID 2
}

// Check a gene is member of the right number of trait modules
TEST_F(PhegeniLoadPheGenISetsFixture, CorrectTraitOwners) {
    EXPECT_EQ(2, modules.entities_to_traits["CFH"].count()) << "The gene is member of the wrong number of traits.";
    EXPECT_TRUE(modules.entities_to_traits["CFH"][phegeni_traits.str_to_int["Macular Degeneration"]]);
    EXPECT_TRUE(modules.entities_to_traits["CFH"][phegeni_traits.str_to_int["Wet Macular Degeneration"]]);

    EXPECT_EQ(1, modules.entities_to_traits["FADS1"].count()) << "The gene is member of the wrong number of traits.";
    EXPECT_TRUE(modules.entities_to_traits["FADS1"][phegeni_traits.str_to_int["Metabolism"]]);
}

class PhegeniConvertModulesWithMapping : public ::testing::Test {

protected:
    virtual void SetUp() {
        gene_modules = loadPheGenIGeneModules(path_file_phegeni, path_file_genes);
        genes = createBimap(path_file_genes);
        proteins = createBimap(path_file_proteins);
        mapping = readMapping(path_file_mapping);
        protein_modules = convertModulesWithMapping(gene_modules,
                                                                  genes,
                                                                  proteins,
                                                                  mapping.second_to_first);
    }

    std::string path_file_phegeni = "../../../../resources/PheGenI/PheGenI_Association_genome_wide_significant_slice.txt";
    std::string path_file_genes = "../../../../resources/Reactome/v70/Genes/genes_slice.csv";
    std::string path_file_proteins = "../../../../resources/Reactome/v70/Proteins/all_proteins_v70.csv";
    std::string path_file_mapping = "../../../../resources/UniProt/mapping_genes_proteins_v70.tab";
    trait_modules gene_modules, protein_modules;
    bimap_str_int genes, proteins;
    entity_mapping mapping;
};

// Return correct trait set names
TEST_F(PhegeniConvertModulesWithMapping, ModuleTraitNamesCorrect) {
    EXPECT_EQ(gene_modules.traits_to_entities.size(), protein_modules.traits_to_entities.size()) << "The gene modules should be the same number of protein modules.";
    EXPECT_TRUE(protein_modules.traits_to_entities.find("Metabolism") != protein_modules.traits_to_entities.end());
    EXPECT_TRUE(protein_modules.traits_to_entities.find("Vascular Endothelial Growth Factors") != protein_modules.traits_to_entities.end());
    EXPECT_TRUE(protein_modules.traits_to_entities.find("Bilirubin") != protein_modules.traits_to_entities.end());
}

// Check module members are correct
TEST_F(PhegeniConvertModulesWithMapping, ModuleTraitSetMembersCorrect){
//    std::cout << "[ ";
//    for(int I = 0; I < proteins.int_to_str.size(); I++){
//        if(protein_modules.traits_to_entities["\"Cholesterol, HDL\""][I]){
//            std::cout << proteins.int_to_str[I] << ", ";
//        }
//    }
//    std::cout << "]\n" << std::endl;
    EXPECT_EQ(2, protein_modules.traits_to_entities["\"Cholesterol, HDL\""].count());
    EXPECT_TRUE(protein_modules.traits_to_entities["\"Cholesterol, HDL\""][proteins.str_to_int["Q15011"]]); // From gene HERPUD1
    EXPECT_TRUE(protein_modules.traits_to_entities["\"Cholesterol, HDL\""][proteins.str_to_int["P11597"]]); // From gene CETP
    EXPECT_FALSE(protein_modules.traits_to_entities["\"Cholesterol, HDL\""][proteins.str_to_int["P02649"]]); // From gene APOE

    EXPECT_EQ(9, protein_modules.traits_to_entities["Bilirubin"].count());
    EXPECT_TRUE(protein_modules.traits_to_entities["Bilirubin"][proteins.str_to_int["P22309"]]); // From gene UGT1A1
    EXPECT_TRUE(protein_modules.traits_to_entities["Bilirubin"][proteins.str_to_int["P22310"]]); // From gene UGT1A4
}

// Return correct protein set names
TEST_F(PhegeniConvertModulesWithMapping, ModuleProteinNamesCorrect) {
    EXPECT_EQ(19, protein_modules.entities_to_traits.size()) << "The proteins in the modules should be NN.";
    EXPECT_TRUE(protein_modules.entities_to_traits.find("P22309") != protein_modules.entities_to_traits.end()); // Comming from gene UGT1A1
    EXPECT_TRUE(protein_modules.entities_to_traits.find("O60427") != protein_modules.entities_to_traits.end()); // Comming from gene FADS1
    EXPECT_TRUE(protein_modules.entities_to_traits.find("Q15011") != protein_modules.entities_to_traits.end());   // Comming from gene HERPUD1
    EXPECT_FALSE(protein_modules.entities_to_traits.find("ALMS1P") != protein_modules.entities_to_traits.end());    // There is no mapping for gene
}

