#include "gtest/gtest.h"
#include <bimap_str_int.hpp>
#include "phegeni.hpp"


class PhegeniLoadPheGenISetsFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        genes = createBimap(path_file_genes);
        auto ret = loadPheGenIGenesAndTraits(path_file_phegeni, genes);
        traits = ret.groups;
        gene_modules = createAndLoadPheGenIGeneModules(path_file_phegeni, genes, traits,
                                                       path_file_gene_interactions,
                                                       "../../../../reports/modules/", ".tsv");
        std::cerr << traits.str_to_int.size() << " === " << traits.int_to_str.size() << "\n";
    }

    std::string path_file_phegeni = "../../Google_tests/resources/PheGenI_Association_genome_wide_significant_slice.txt";
    std::string path_file_genes = "../../Google_tests/resources/genes_slice.csv";
    std::string path_file_gene_interactions = "../../Google_tests/resources/gene_interactions.tab";
    bimap_str_int genes, traits;
    modules gene_modules;
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
    EXPECT_EQ(8, traits.int_to_str.size());
    EXPECT_EQ(8, traits.str_to_int.size());

    // Check that the genes are correct
    EXPECT_EQ(25, genes.int_to_str.size());
    EXPECT_EQ(25, genes.str_to_int.size());

    // Check genes are sorted
    EXPECT_EQ(1, traits.str_to_int.at("\"Cholesterol, LDL\""));
    EXPECT_EQ("\"Cholesterol, LDL\"", traits.int_to_str[1]);
    EXPECT_EQ(2, traits.str_to_int.at("Bilirubin"));
    EXPECT_EQ("Bilirubin", traits.int_to_str[2]);
    EXPECT_EQ(7, traits.str_to_int.at("Wet Macular Degeneration"));
    EXPECT_EQ("Wet Macular Degeneration", traits.int_to_str[7]);
}

// Test the trait to gene bitsets contains the right number of trait keys
// Test the gene to trait bitsets contains the right number of gene keys
TEST_F(PhegeniLoadPheGenISetsFixture, TraitKeys) {
    EXPECT_EQ(8, gene_modules.group_to_members.size()) << "It has the wrong number of Traits as keys of the map.";
    EXPECT_EQ(25, gene_modules.member_to_groups.size()) << "It has the wrong number of GENES as keys of the map.";
}

// Test the trait to gene bitset contains the correct number of bits in the bitsets
// Test the gene to trait bitsets contain the correct number of bits in the bitsets
TEST_F(PhegeniLoadPheGenISetsFixture, BitsetSizes) {
    EXPECT_EQ(25, gene_modules.group_to_members["Bilirubin"].size()) << "The bitsets size are not the number of genes.";
    EXPECT_EQ(8, gene_modules.member_to_groups["UGT1A1"].size()) << "The bitsets size are not the number of traits.";
}

// Check a trait has the right number of genes as members
TEST_F(PhegeniLoadPheGenISetsFixture, CorrectGeneMembers) {
    EXPECT_EQ(9, gene_modules.group_to_members["Bilirubin"].count());
    EXPECT_TRUE(gene_modules.group_to_members["Bilirubin"][genes.str_to_int["UGT1A1"]]); // From column GENE ID
    EXPECT_TRUE(gene_modules.group_to_members["Bilirubin"][genes.str_to_int["UGT1A4"]]); // From column GENE ID 2
    EXPECT_FALSE(gene_modules.group_to_members["Bilirubin"][genes.str_to_int["AAAA"]]);  // From column GENE ID
    EXPECT_FALSE(gene_modules.group_to_members["Bilirubin"][genes.str_to_int["DDDD"]]);  // From column GENE ID 2

    EXPECT_EQ(2, gene_modules.group_to_members["\"Cholesterol, HDL\""].count());
    EXPECT_TRUE(
            gene_modules.group_to_members["\"Cholesterol, HDL\""][genes.str_to_int["HERPUD1"]]); // From column GENE ID
    EXPECT_TRUE(
            gene_modules.group_to_members["\"Cholesterol, HDL\""][genes.str_to_int["CETP"]]); // From column GENE ID 2
    EXPECT_FALSE(gene_modules.group_to_members["\"Cholesterol, HDL\""][genes.str_to_int["APOE"]]);

    EXPECT_EQ(9, gene_modules.group_to_members["Bilirubin"].count());
    EXPECT_TRUE(gene_modules.group_to_members["Bilirubin"][genes.str_to_int["UGT1A1"]]); // From column GENE ID
    EXPECT_TRUE(gene_modules.group_to_members["Bilirubin"][genes.str_to_int["UGT1A4"]]); // From column GENE ID 2
    EXPECT_FALSE(gene_modules.group_to_members["Bilirubin"][genes.str_to_int["AAAA"]]);  // From column GENE ID
    EXPECT_FALSE(gene_modules.group_to_members["Bilirubin"][genes.str_to_int["DDDD"]]);  // From column GENE ID 2
}

// Check a gene is member of the right number of trait modules
TEST_F(PhegeniLoadPheGenISetsFixture, CorrectTraitOwners) {
    std::cout << "[ ";
    for (int I = 0; I < gene_modules.member_to_groups["CFH"].size(); I++) {
        if (gene_modules.member_to_groups["CFH"][I]) {
            std::cout << traits.int_to_str[I] << ", ";
        }
    }
    std::cout << "]\n" << std::endl;

    EXPECT_EQ(2, gene_modules.member_to_groups["CFH"].count()) << "The gene is member of the wrong number of traits.";
    EXPECT_TRUE(gene_modules.member_to_groups["CFH"][traits.str_to_int["Macular Degeneration"]]);
    EXPECT_TRUE(gene_modules.member_to_groups["CFH"][traits.str_to_int["Wet Macular Degeneration"]]);

    EXPECT_EQ(1, gene_modules.member_to_groups["FADS1"].count()) << "The gene is member of the wrong number of traits.";
    EXPECT_TRUE(gene_modules.member_to_groups["FADS1"][traits.str_to_int["Metabolism"]]);
}

class PhegeniConvertModulesWithMapping : public ::testing::Test {

protected:
    virtual void SetUp() {
        genes = createBimap(path_file_genes);
        proteins = createBimap(path_file_proteins);
        mapping = readMapping(path_file_mapping);
        auto ret = loadPheGenIGenesAndTraits(path_file_phegeni, genes);
        traits = ret.groups;
        gene_modules = createAndLoadPheGenIGeneModules(path_file_phegeni, genes, traits, path_file_gene_interactions,
                                                       "../../../../reports/modules/", ".tsv");
        protein_modules = convertModulesWithMapping(gene_modules,
                                                    genes,
                                                    proteins,
                                                    traits,
                                                    mapping.second_to_first);
    }

    std::string path_file_phegeni = "../../Google_tests/resources/PheGenI_Association_genome_wide_significant_slice.txt";
    std::string path_file_genes = "../../Google_tests/resources/genes_slice.csv";
    std::string path_file_gene_interactions = "../../Google_tests/resources/gene_interactions.tab";
    std::string path_file_proteins = "../../Google_tests/resources/proteins_slice.csv";
    std::string path_file_mapping = "../../../../resources/UniProt/mapping_proteins_to_genes_v70.tab";
    modules gene_modules, protein_modules;
    bimap_str_int genes, proteins, phegeni_genes, traits;
    bidirectional_mapping mapping;
};

// Return correct trait set names
TEST_F(PhegeniConvertModulesWithMapping, ModuleTraitNamesCorrect) {
    EXPECT_EQ(gene_modules.group_to_members.size(), protein_modules.group_to_members.size())
                        << "The gene modules should be the same number of protein modules.";
    EXPECT_TRUE(protein_modules.group_to_members.find("Metabolism") != protein_modules.group_to_members.end());
    EXPECT_TRUE(protein_modules.group_to_members.find("Vascular Endothelial Growth Factors") !=
                protein_modules.group_to_members.end());
    EXPECT_TRUE(protein_modules.group_to_members.find("Bilirubin") != protein_modules.group_to_members.end());
}

// Check module members are correct
TEST_F(PhegeniConvertModulesWithMapping, ModuleTraitSetMembersCorrect) {
    EXPECT_EQ(2, protein_modules.group_to_members["\"Cholesterol, HDL\""].count());
    EXPECT_TRUE(
            protein_modules.group_to_members["\"Cholesterol, HDL\""][proteins.str_to_int["Q15011"]]); // From gene HERPUD1
    EXPECT_TRUE(
            protein_modules.group_to_members["\"Cholesterol, HDL\""][proteins.str_to_int["P11597"]]); // From gene CETP
    EXPECT_FALSE(
            protein_modules.group_to_members["\"Cholesterol, HDL\""][proteins.str_to_int["P02649"]]); // From gene APOE

    EXPECT_EQ(9, protein_modules.group_to_members["Bilirubin"].count());
    EXPECT_TRUE(protein_modules.group_to_members["Bilirubin"][proteins.str_to_int["P22309"]]); // From gene UGT1A1
    EXPECT_TRUE(protein_modules.group_to_members["Bilirubin"][proteins.str_to_int["P22310"]]); // From gene UGT1A4
}

// Return correct protein set names
TEST_F(PhegeniConvertModulesWithMapping, ModuleProteinNamesCorrect) {
    EXPECT_EQ(19, protein_modules.member_to_groups.size()) << "The proteins in the modules should be 19.";
    EXPECT_TRUE(protein_modules.member_to_groups.find("P22309") !=
                protein_modules.member_to_groups.end()); // Comming from gene UGT1A1
    EXPECT_TRUE(protein_modules.member_to_groups.find("O60427") !=
                protein_modules.member_to_groups.end()); // Comming from gene FADS1
    EXPECT_TRUE(protein_modules.member_to_groups.find("Q15011") !=
                protein_modules.member_to_groups.end());   // Comming from gene HERPUD1
    EXPECT_FALSE(protein_modules.member_to_groups.find("ALMS1P") !=
                 protein_modules.member_to_groups.end());    // There is no mapping for gene
}

// Check protein owners are correct
TEST_F(PhegeniConvertModulesWithMapping, ModuleProteinOwnerSetCorrect) {
    EXPECT_EQ(traits.int_to_str.size(), protein_modules.member_to_groups.begin()->second.size())
                        << "The bitset size should be the number of traits.";

    EXPECT_EQ(1, protein_modules.member_to_groups["P11597"].count());
    EXPECT_TRUE(protein_modules.member_to_groups["Q15011"][traits.str_to_int["\"Cholesterol, HDL\""]]);

    EXPECT_EQ(1, protein_modules.member_to_groups["P22309"].count());
    EXPECT_TRUE(protein_modules.member_to_groups["P22309"][traits.str_to_int["Bilirubin"]]);
}

class CreatePheGenIModulesFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        genes = createBimap(path_file_genes);
        proteins = createBimap(path_file_proteins);
        proteoforms = createBimap(path_file_proteoforms);
        auto ret = loadPheGenIGenesAndTraits(path_file_phegeni_modules, genes);
        traits = ret.groups;

        gene_modules = createAndLoadPheGenIGeneModules(path_file_phegeni_modules, genes, traits,
                                                       path_file_gene_interactions, "../../../../reports/modules/",
                                                       ".tsv");

        bidirectional_mapping mapping_genes_to_proteins = readMapping(path_file_mapping_genes_to_proteins,
                                                                      true, false);
        protein_modules = createAndLoadPheGenIModules(gene_modules, genes, proteins, traits,
                                                      mapping_genes_to_proteins.first_to_second,
                                                      path_file_protein_interactions,
                                                      "../../../../reports/modules/", "proteins", ".tsv");

        bidirectional_mapping mapping_proteins_to_proteoforms = readMapping(path_file_mapping_proteins_to_proteoforms,
                                                                            true,
                                                                            true);
        proteoform_modules = createAndLoadPheGenIModules(protein_modules, proteins, proteoforms, traits,
                                                         mapping_proteins_to_proteoforms.first_to_second,
                                                         path_file_proteoform_interactions,
                                                         "../../../../reports/modules/", "proteoforms", ".tsv");
    }

    bimap_str_int genes, proteins, proteoforms, traits;
    modules gene_modules, protein_modules, proteoform_modules;

    std::string path_file_phegeni_modules = "../../Google_tests/resources/createPheGenIModules/modules.csv";
    std::string path_file_genes = "../../Google_tests/resources/createPheGenIModules/genes.csv";
    std::string path_file_proteins = "../../Google_tests/resources/createPheGenIModules/proteins.csv";
    std::string path_file_proteoforms = "../../Google_tests/resources/createPheGenIModules/proteoforms.csv";
    std::string path_file_gene_interactions = "../../Google_tests/resources/createPheGenIModules/gene_interactions.csv";
    std::string path_file_protein_interactions = "../../Google_tests/resources/createPheGenIModules/protein_interactions.csv";
    std::string path_file_proteoform_interactions = "../../Google_tests/resources/createPheGenIModules/proteoform_interactions.csv";
    std::string path_file_mapping_genes_to_proteins = "../../Google_tests/resources/createPheGenIModules/mapping_genes_to_proteins.csv";
    std::string path_file_mapping_proteins_to_proteoforms = "../../Google_tests/resources/createPheGenIModules/mapping_proteins_to_proteoforms.csv";
};

TEST_F(CreatePheGenIModulesFixture, CorrectTraits) {
    EXPECT_EQ(2, traits.str_to_int.size());
    EXPECT_EQ(2, traits.int_to_str.size());

    EXPECT_EQ(0, traits.str_to_int["GROUP1"]);
    EXPECT_EQ(1, traits.str_to_int["GROUP2"]);

    EXPECT_EQ("GROUP1", traits.int_to_str[0]);
    EXPECT_EQ("GROUP2", traits.int_to_str[1]);
}

// Check genes, proteins and proteoform bimaps are correct
TEST_F(CreatePheGenIModulesFixture, CorrectGenes) {
    EXPECT_EQ(10, genes.int_to_str.size());
    EXPECT_EQ(10, genes.str_to_int.size());

    EXPECT_FALSE(hasKey(genes.str_to_int, (std::string) "GENES"));
    EXPECT_FALSE(hasValue(genes.int_to_str, (std::string) "GENES"));

    EXPECT_EQ("A", genes.int_to_str[0]);
    EXPECT_EQ("C", genes.int_to_str[2]);
    EXPECT_EQ("G", genes.int_to_str[6]);
    EXPECT_EQ("J", genes.int_to_str[9]);

    EXPECT_EQ(0, genes.str_to_int["A"]);
    EXPECT_EQ(1, genes.str_to_int["B"]);
    EXPECT_EQ(3, genes.str_to_int["D"]);
    EXPECT_EQ(4, genes.str_to_int["E"]);
    EXPECT_EQ(8, genes.str_to_int["I"]);
}

TEST_F(CreatePheGenIModulesFixture, CorrectProteins) {
    EXPECT_EQ(20, proteins.int_to_str.size());
    EXPECT_EQ(20, proteins.str_to_int.size());

    EXPECT_FALSE(hasKey(proteins.str_to_int, (std::string) "PROTEINS"));
    EXPECT_FALSE(hasValue(proteins.int_to_str, (std::string) "PROTEINS"));

    EXPECT_EQ("A1", proteins.int_to_str[0]);
    EXPECT_EQ("B1", proteins.int_to_str[3]);
    EXPECT_EQ("D2", proteins.int_to_str[8]);
    EXPECT_EQ("I1", proteins.int_to_str[16]);
    EXPECT_EQ("J3", proteins.int_to_str[19]);

    EXPECT_EQ(0, proteins.str_to_int["A1"]);
    EXPECT_EQ(4, proteins.str_to_int["C1"]);
    EXPECT_EQ(9, proteins.str_to_int["E1"]);
    EXPECT_EQ(19, proteins.str_to_int["J3"]);
}

TEST_F(CreatePheGenIModulesFixture, CorrectProteoforms) {
    EXPECT_EQ(26, proteoforms.int_to_str.size());
    EXPECT_EQ(26, proteoforms.str_to_int.size());

    EXPECT_FALSE(hasKey(proteoforms.str_to_int, (std::string) "PROTEOFORMS"));
    EXPECT_FALSE(hasValue(proteoforms.int_to_str, (std::string) "PROTEOFORMS"));

    EXPECT_EQ("A1_1", proteoforms.int_to_str[0]);
    EXPECT_EQ("C1_1", proteoforms.int_to_str[5]);
    EXPECT_EQ("D1_2", proteoforms.int_to_str[10]);
    EXPECT_EQ("F2_1", proteoforms.int_to_str[15]);
    EXPECT_EQ("J3_1", proteoforms.int_to_str[25]);

    EXPECT_EQ(0, proteoforms.str_to_int["A1_1"]);
    EXPECT_EQ(2, proteoforms.str_to_int["A1_3"]);
    EXPECT_EQ(7, proteoforms.str_to_int["C2_1"]);
    EXPECT_EQ(12, proteoforms.str_to_int["E1_1"]);
    EXPECT_EQ(19, proteoforms.str_to_int["G2_1"]);
    EXPECT_EQ(25, proteoforms.str_to_int["J3_1"]);
}

// Check members of modules are correct

// Check gene modules are correct
TEST_F(CreatePheGenIModulesFixture, CorrectGeneModules) {
    // Check there are 2 modules
    EXPECT_EQ(2, getKeys(gene_modules.group_to_members).size());
    EXPECT_EQ(2, gene_modules.member_to_groups.begin()->second.size());

    // Check there are 10 genes
    EXPECT_EQ(10, getKeys(gene_modules.member_to_groups).size());
    EXPECT_EQ(10, gene_modules.group_to_members.begin()->second.size());

    // Check correct members
    EXPECT_EQ(3, gene_modules.group_to_members.at("GROUP1").count());
    EXPECT_TRUE(gene_modules.group_to_members.at("GROUP1")[genes.str_to_int["A"]]);
    EXPECT_TRUE(gene_modules.group_to_members.at("GROUP1")[genes.str_to_int["C"]]);

    EXPECT_EQ(4, gene_modules.group_to_members.at("GROUP2").count());
    EXPECT_TRUE(gene_modules.group_to_members.at("GROUP2")[genes.str_to_int["D"]]);
    EXPECT_TRUE(gene_modules.group_to_members.at("GROUP2")[genes.str_to_int["E"]]);
    EXPECT_TRUE(gene_modules.group_to_members.at("GROUP2")[genes.str_to_int["F"]]);
    EXPECT_TRUE(gene_modules.group_to_members.at("GROUP2")[genes.str_to_int["G"]]);

    // Check correct owners
    EXPECT_EQ(1, gene_modules.member_to_groups.at("A").count());
    EXPECT_EQ(1, gene_modules.member_to_groups.at("C").count());
    EXPECT_EQ(1, gene_modules.member_to_groups.at("E").count());
    EXPECT_EQ(1, gene_modules.member_to_groups.at("G").count());
    EXPECT_EQ(0, gene_modules.member_to_groups.at("H").count());
    EXPECT_EQ(0, gene_modules.member_to_groups.at("J").count());

    EXPECT_TRUE(gene_modules.member_to_groups.at("B")[traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(gene_modules.member_to_groups.at("F")[traits.str_to_int["GROUP2"]]);
    EXPECT_FALSE(gene_modules.member_to_groups.at("I")[traits.str_to_int["GROUP1"]]);
}

// Check protein modules are correct
TEST_F(CreatePheGenIModulesFixture, CorrectProteinModulesSizes) {

    // -- Check the number of modules is the same
    // Check there are still 2 modules
    EXPECT_EQ(2, getKeys(protein_modules.group_to_members).size());
    EXPECT_EQ(2, protein_modules.member_to_groups.begin()->second.size());

    // Check there are 20 proteins
    EXPECT_EQ(20, getKeys(protein_modules.member_to_groups).size());
    EXPECT_EQ(20, protein_modules.group_to_members.begin()->second.size());
}

TEST_F(CreatePheGenIModulesFixture, CorrectProteinModulesMembers) {
    // Check correct members
    EXPECT_EQ(4, protein_modules.group_to_members["GROUP1"].count());
    // -- Check a gene with two protein products appears in the new modules
    EXPECT_TRUE(protein_modules.group_to_members["GROUP1"][proteins.str_to_int["A1"]]);
    // -- Check that one of the protein products of a gene gets removed of the modules, because that one does not interact with other members
    EXPECT_FALSE(protein_modules.group_to_members["GROUP1"][proteins.str_to_int["A2"]]);
    // -- Check that one of the protein products of a gene gets removed of the modules, because that one does not interact with other members
    EXPECT_FALSE(protein_modules.group_to_members["GROUP1"][proteins.str_to_int["A3"]]);
    EXPECT_TRUE(protein_modules.group_to_members["GROUP1"][proteins.str_to_int["B1"]]);
    // -- Check a gene with two protein products appears in the new modules
    EXPECT_TRUE(protein_modules.group_to_members["GROUP1"][proteins.str_to_int["C1"]]);
    // -- Check that one of the protein products of a gene gets removed of the modules, because that one does not interact with other members
    EXPECT_FALSE(protein_modules.group_to_members["GROUP1"][proteins.str_to_int["C2"]]);
    EXPECT_TRUE(protein_modules.group_to_members["GROUP1"][proteins.str_to_int["C3"]]);
    EXPECT_FALSE(protein_modules.group_to_members["GROUP1"][proteins.str_to_int["F1"]]);

    EXPECT_EQ(6, protein_modules.group_to_members["GROUP2"].count());
    EXPECT_TRUE(protein_modules.group_to_members["GROUP2"][proteins.str_to_int["D1"]]);
    EXPECT_TRUE(protein_modules.group_to_members["GROUP2"][proteins.str_to_int["D2"]]);
    EXPECT_TRUE(protein_modules.group_to_members["GROUP2"][proteins.str_to_int["E1"]]);
    EXPECT_TRUE(protein_modules.group_to_members["GROUP2"][proteins.str_to_int["F1"]]);
    EXPECT_FALSE(protein_modules.group_to_members["GROUP2"][proteins.str_to_int["F2"]]);
    EXPECT_TRUE(protein_modules.group_to_members["GROUP2"][proteins.str_to_int["G1"]]);
    EXPECT_TRUE(protein_modules.group_to_members["GROUP2"][proteins.str_to_int["G2"]]);
}

TEST_F(CreatePheGenIModulesFixture, CorrectProteinModulesOwners) {
    // -- -- Check the proteins have/do not have ownership of the traits
    EXPECT_EQ(1, protein_modules.member_to_groups["A1"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["A2"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["A3"].count());
    EXPECT_EQ(1, protein_modules.member_to_groups["B1"].count());
    EXPECT_EQ(1, protein_modules.member_to_groups["C1"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["C2"].count());
    EXPECT_EQ(1, protein_modules.member_to_groups["C3"].count());

    EXPECT_TRUE(protein_modules.member_to_groups["A1"][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(protein_modules.member_to_groups["B1"][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(protein_modules.member_to_groups["C1"][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(protein_modules.member_to_groups["C3"][traits.str_to_int["GROUP1"]]);

    EXPECT_EQ(1, protein_modules.member_to_groups["D1"].count());
    EXPECT_EQ(1, protein_modules.member_to_groups["D2"].count());
    EXPECT_EQ(1, protein_modules.member_to_groups["E1"].count());
    EXPECT_EQ(1, protein_modules.member_to_groups["F1"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["F2"].count());
    EXPECT_EQ(1, protein_modules.member_to_groups["G1"].count());
    EXPECT_EQ(1, protein_modules.member_to_groups["G2"].count());

    EXPECT_TRUE(protein_modules.member_to_groups["D1"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups["D2"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups["E1"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups["F1"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups["G1"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups["G2"][traits.str_to_int["GROUP2"]]);

    EXPECT_EQ(0, protein_modules.member_to_groups["H1"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["H2"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["I1"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["J1"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["J2"].count());
    EXPECT_EQ(0, protein_modules.member_to_groups["J3"].count());
}

// Check proteoform modules are correct
// -- Check the number of modules is the same
// -- Check a gene with two protein products appears in the new modules
// -- -- Check the proteins have/do not have ownership of the traits
// -- Check that one of the protein products of a gene gets removed of the modules, because that one does not interact with other members
// -- -- Check the proteins have/do not have ownership of the traits
// -- Check random members are correct
// -- -- Check the inverse ownership is correct
TEST_F(CreatePheGenIModulesFixture, CorrectProteoformModulesSizes) {

    // -- Check the number of modules is the same
    // Check there are still 2 modules
    EXPECT_EQ(2, getKeys(proteoform_modules.group_to_members).size());
    EXPECT_EQ(2, proteoform_modules.member_to_groups.begin()->second.size());

    // Check there are 20 proteins
    EXPECT_EQ(26, getKeys(proteoform_modules.member_to_groups).size());
    EXPECT_EQ(26, proteoform_modules.group_to_members.begin()->second.size());
}

TEST_F(CreatePheGenIModulesFixture, CorrectProteoformModulesMembers) {

    EXPECT_EQ(5, proteoform_modules.group_to_members["GROUP1"].count());
    EXPECT_FALSE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["A1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["A1_2"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["A1_3"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["B1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["B1_2"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["C1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["C1_2"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["C2_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP1"][proteoforms.str_to_int["C3_1"]]);

    EXPECT_EQ(8, proteoform_modules.group_to_members["GROUP2"].count());
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["D1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["D1_2"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["D2_1"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["E1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["E1_2"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["F1_1"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["F2_1"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["F2_2"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["G1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["G1_2"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members["GROUP2"][proteoforms.str_to_int["G2_1"]]);
}

TEST_F(CreatePheGenIModulesFixture, CorrectProtoformModulesOwners) {
    EXPECT_EQ(0, proteoform_modules.member_to_groups["A1_1"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["A1_2"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["A1_3"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["B1_1"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["B1_2"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["C1_1"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["C1_2"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["C2_1"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["C3_1"].count());

    EXPECT_TRUE(proteoform_modules.member_to_groups["A1_2"][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["B1_1"][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["B1_2"][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["C1_2"][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["C3_1"][traits.str_to_int["GROUP1"]]);

    EXPECT_EQ(1, proteoform_modules.member_to_groups["D1_1"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["D1_2"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["D2_1"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["E1_1"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["E1_2"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["F1_1"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["F2_1"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["F2_2"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["G1_1"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["G1_2"].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups["G2_1"].count());

    EXPECT_TRUE(proteoform_modules.member_to_groups["D1_1"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["D1_2"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["D2_1"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["E1_2"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["F1_1"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["G1_1"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["G1_2"][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups["G2_1"][traits.str_to_int["GROUP2"]]);

    EXPECT_EQ(0, proteoform_modules.member_to_groups["H1_1"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["H2_1"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["I1_1"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["J1_1"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["J2_1"].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups["J3_1"].count());
}

