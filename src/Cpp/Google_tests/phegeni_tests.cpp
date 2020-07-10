#include "gtest/gtest.h"
#include <bimap_str_int.hpp>
#include "phegeni.hpp"


class PhegeniLoadPheGenISetsFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        genes = createBimap(path_file_genes);
        auto ret = loadPheGenIGenesAndTraits(path_file_phegeni, genes);
        traits = ret.groups;
        gene_modules = createGeneModules(path_file_phegeni, genes, traits,
                                         path_file_gene_interactions,
                                         "../../../../reports/modules/", ".tsv", false);
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
    int trait_index = traits.str_to_int["Bilirubin"];
    EXPECT_EQ(25, gene_modules.group_to_members[trait_index].size()) << "The bitsets size are not the number of genes.";

    int gene_index = genes.str_to_int["UGT1A1"];
    EXPECT_EQ(8, gene_modules.member_to_groups[gene_index].size()) << "The bitsets size are not the number of traits.";
}

// Check a trait has the right number of genes as members
TEST_F(PhegeniLoadPheGenISetsFixture, CorrectGeneMembers) {
    int gene_index, trait_index;

    trait_index = traits.str_to_int["Bilirubin"];
    EXPECT_EQ(9, gene_modules.group_to_members[trait_index].count());
    EXPECT_TRUE(gene_modules.group_to_members[trait_index][genes.str_to_int["UGT1A1"]]); // From column GENE ID
    EXPECT_TRUE(gene_modules.group_to_members[trait_index][genes.str_to_int["UGT1A4"]]); // From column GENE ID 2
    EXPECT_FALSE(gene_modules.group_to_members[trait_index][genes.str_to_int["AAAA"]]);  // From column GENE ID
    EXPECT_FALSE(gene_modules.group_to_members[trait_index][genes.str_to_int["DDDD"]]);  // From column GENE ID 2

    trait_index = traits.str_to_int["\"Cholesterol, HDL\""];
    EXPECT_EQ(2, gene_modules.group_to_members[trait_index].count());
    EXPECT_TRUE(gene_modules.group_to_members[trait_index][genes.str_to_int["HERPUD1"]]); // From column GENE ID
    EXPECT_TRUE(gene_modules.group_to_members[trait_index][genes.str_to_int["CETP"]]); // From column GENE ID 2
    EXPECT_FALSE(gene_modules.group_to_members[trait_index][genes.str_to_int["APOE"]]);

    trait_index = traits.str_to_int["Bilirubin"];
    EXPECT_EQ(9, gene_modules.group_to_members[trait_index].count());
    EXPECT_TRUE(gene_modules.group_to_members[trait_index][genes.str_to_int["UGT1A1"]]); // From column GENE ID
    EXPECT_TRUE(gene_modules.group_to_members[trait_index][genes.str_to_int["UGT1A4"]]); // From column GENE ID 2
    EXPECT_FALSE(gene_modules.group_to_members[trait_index][genes.str_to_int["AAAA"]]);  // From column GENE ID
    EXPECT_FALSE(gene_modules.group_to_members[trait_index][genes.str_to_int["DDDD"]]);  // From column GENE ID 2
}

// Check a gene is member of the right number of trait modules
TEST_F(PhegeniLoadPheGenISetsFixture, CorrectTraitOwners) {
    int gene_index = genes.str_to_int["CFH"];
    std::cout << "[ ";
    for (int trait_index = 0; trait_index < gene_modules.member_to_groups[gene_index].size(); trait_index++) {
        if (gene_modules.member_to_groups[gene_index][trait_index]) {
            std::cout << traits.int_to_str[trait_index] << ", ";
        }
    }
    std::cout << "]\n" << std::endl;

    gene_index = genes.str_to_int["CFH"];
    EXPECT_EQ(2, gene_modules.member_to_groups[gene_index].count())
                        << "The gene is member of the wrong number of traits.";
    EXPECT_TRUE(gene_modules.member_to_groups[gene_index][traits.str_to_int["Macular Degeneration"]]);
    EXPECT_TRUE(gene_modules.member_to_groups[gene_index][traits.str_to_int["Wet Macular Degeneration"]]);

    gene_index = genes.str_to_int["FADS1"];
    EXPECT_EQ(1, gene_modules.member_to_groups[gene_index].count())
                        << "The gene is member of the wrong number of traits.";
    EXPECT_TRUE(gene_modules.member_to_groups[gene_index][traits.str_to_int["Metabolism"]]);
}

class PhegeniConvertModulesWithMapping : public ::testing::Test {

protected:
    virtual void SetUp() {
        genes = createBimap(path_file_genes);
        proteins = createBimap(path_file_proteins);
        mapping = readMapping(path_file_mapping);
        auto ret = loadPheGenIGenesAndTraits(path_file_phegeni, genes);
        traits = ret.groups;
        gene_modules = createGeneModules(path_file_phegeni, genes, traits, path_file_gene_interactions,
                                         "../../../../reports/modules/", ".tsv", false);
        protein_modules = convertModulesWithMapping(gene_modules,
                                                    genes,
                                                    proteins,
                                                    traits,
                                                    mapping.second_to_first);
        writeModulesSingleFile("../../../../reports/modules/", "proteins", ".tsv", protein_modules, traits, proteins);
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
    EXPECT_EQ(protein_modules.group_to_members.size(), traits.int_to_str.size());
}

// Check module members are correct
TEST_F(PhegeniConvertModulesWithMapping, ModuleTraitSetMembersCorrect) {
    int trait_index = traits.str_to_int["\"Cholesterol, HDL\""];
    EXPECT_EQ(2, protein_modules.group_to_members[trait_index].count());
    EXPECT_TRUE(
            protein_modules.group_to_members[trait_index][proteins.str_to_int["Q15011"]]); // From gene HERPUD1
    EXPECT_TRUE(
            protein_modules.group_to_members[trait_index][proteins.str_to_int["P11597"]]); // From gene CETP
    EXPECT_FALSE(
            protein_modules.group_to_members[trait_index][proteins.str_to_int["P02649"]]); // From gene APOE

    trait_index = traits.str_to_int["Bilirubin"];
    std::cout << "There are " << gene_modules.group_to_members.size() << " gene modules.\n";
    EXPECT_EQ(9, gene_modules.group_to_members[trait_index].count());
    for (int member_index = 0; member_index < genes.str_to_int.size(); member_index++) {
        if (gene_modules.group_to_members[trait_index][member_index])
            std::cout << genes.int_to_str[member_index] << ", ";
    }
    std::cout << std::endl;
    std::cout << "There are " << gene_modules.group_to_members.size() << " protein modules.\n";
    EXPECT_EQ(9, protein_modules.group_to_members[trait_index].count());
    for (int member_index = 0; member_index < proteins.str_to_int.size(); member_index++) {
        if (protein_modules.group_to_members[trait_index][member_index])
            std::cout << genes.int_to_str[member_index] << ", ";
    }
    std::cout << std::endl;
    EXPECT_TRUE(protein_modules.group_to_members[trait_index][proteins.str_to_int["P22309"]]); // From gene UGT1A1
    EXPECT_TRUE(protein_modules.group_to_members[trait_index][proteins.str_to_int["P22310"]]); // From gene UGT1A4
}

// Return correct protein set names
TEST_F(PhegeniConvertModulesWithMapping, ModuleProteinNamesCorrect) {
    EXPECT_EQ(19, protein_modules.member_to_groups.size()) << "The proteins in the modules should be 19.";
    EXPECT_TRUE(hasKey(proteins.str_to_int, static_cast<std::string>("P22309"))); // Comming from gene UGT1A1
    EXPECT_TRUE(hasKey(proteins.str_to_int, static_cast<std::string>("Q15011")));   // Comming from gene FADS1
    EXPECT_FALSE(hasKey(proteins.str_to_int, static_cast<std::string>("ALMS1P")));    // There is no mapping for gene
}

// Check protein owners are correct
TEST_F(PhegeniConvertModulesWithMapping, ModuleProteinOwnerSetCorrect) {
    EXPECT_EQ(traits.int_to_str.size(), protein_modules.member_to_groups[0].size())
                        << "The bitset size should be the number of traits.";

    int protein_index = proteins.str_to_int["P11597"];
    EXPECT_EQ(1, protein_modules.member_to_groups[protein_index].count());

    protein_index = proteins.str_to_int["Q15011"];
    int trait_index = traits.str_to_int["\"Cholesterol, HDL\""];
    EXPECT_TRUE(protein_modules.member_to_groups[protein_index][trait_index]);

    protein_index = proteins.str_to_int["P22309"];
    trait_index = traits.str_to_int["Bilirubin"];
    EXPECT_EQ(1, protein_modules.member_to_groups[protein_index].count());
    EXPECT_TRUE(protein_modules.member_to_groups[protein_index][trait_index]);
}

class CreatePheGenIModulesFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        genes = createBimap(path_file_genes);
        proteins = createBimap(path_file_proteins);
        proteoforms = createBimap(path_file_proteoforms);
        auto ret = loadPheGenIGenesAndTraits(path_file_phegeni_modules, genes);
        traits = ret.groups;

        gene_modules = createGeneModules(path_file_phegeni_modules, genes, traits,
                                         path_file_gene_interactions, "../../../../reports/modules/",
                                         ".tsv", false);

        bidirectional_mapping mapping_genes_to_proteins = readMapping(path_file_mapping_genes_to_proteins,
                                                                      true, false);
        protein_modules = createProteinOrProteoformModules(gene_modules, genes, proteins, traits,
                                                           mapping_genes_to_proteins.first_to_second,
                                                           path_file_protein_interactions,
                                                           "../../../../reports/modules/", "proteins", ".tsv",
                                                           false);

        bidirectional_mapping mapping_proteins_to_proteoforms = readMapping(path_file_mapping_proteins_to_proteoforms,
                                                                            true,
                                                                            true);
        proteoform_modules = createProteinOrProteoformModules(protein_modules, proteins, proteoforms, traits,
                                                              mapping_proteins_to_proteoforms.first_to_second,
                                                              path_file_proteoform_interactions,
                                                              "../../../../reports/modules/", "proteoforms",
                                                              ".tsv", false);
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
    EXPECT_EQ(2, gene_modules.group_to_members.size());
    EXPECT_EQ(2, gene_modules.member_to_groups[0].size());

    // Check there are 10 genes
    EXPECT_EQ(10, gene_modules.member_to_groups.size());
    EXPECT_EQ(10, gene_modules.group_to_members[0].size());

    // Check correct members
    // GROUP1
    EXPECT_EQ(3, gene_modules.group_to_members[0].count());
    EXPECT_TRUE(gene_modules.group_to_members[0][genes.str_to_int["A"]]);
    EXPECT_TRUE(gene_modules.group_to_members[0][genes.str_to_int["C"]]);

    // GROUP2
    EXPECT_EQ(4, gene_modules.group_to_members[1].count());
    EXPECT_TRUE(gene_modules.group_to_members[1][genes.str_to_int["D"]]);
    EXPECT_TRUE(gene_modules.group_to_members[1][genes.str_to_int["E"]]);
    EXPECT_TRUE(gene_modules.group_to_members[1][genes.str_to_int["F"]]);
    EXPECT_TRUE(gene_modules.group_to_members[1][genes.str_to_int["G"]]);

    // Check correct owners
    int gene_index = genes.str_to_int["A"];
    EXPECT_EQ(1, gene_modules.member_to_groups[gene_index].count());
    gene_index = genes.str_to_int["C"];
    EXPECT_EQ(1, gene_modules.member_to_groups[gene_index].count());
    gene_index = genes.str_to_int["E"];
    EXPECT_EQ(1, gene_modules.member_to_groups[gene_index].count());
    gene_index = genes.str_to_int["G"];
    EXPECT_EQ(1, gene_modules.member_to_groups[gene_index].count());
    gene_index = genes.str_to_int["H"];
    EXPECT_EQ(0, gene_modules.member_to_groups[gene_index].count());
    gene_index = genes.str_to_int["J"];
    EXPECT_EQ(0, gene_modules.member_to_groups[gene_index].count());

    gene_index = genes.str_to_int["B"];
    EXPECT_TRUE(gene_modules.member_to_groups[gene_index][traits.str_to_int["GROUP1"]]);
    gene_index = genes.str_to_int["F"];
    EXPECT_TRUE(gene_modules.member_to_groups[gene_index][traits.str_to_int["GROUP2"]]);
    gene_index = genes.str_to_int["I"];
    EXPECT_FALSE(gene_modules.member_to_groups[gene_index][traits.str_to_int["GROUP1"]]);
}

// Check protein modules are correct
TEST_F(CreatePheGenIModulesFixture, CorrectProteinModulesSizes) {

    // -- Check the number of modules is the same
    // Check there are still 2 modules
    EXPECT_EQ(2, protein_modules.group_to_members.size());
    EXPECT_EQ(2, protein_modules.member_to_groups[0].size());

    // Check there are 20 proteins
    EXPECT_EQ(20, protein_modules.member_to_groups.size());
    EXPECT_EQ(20, protein_modules.group_to_members[0].size());
}

TEST_F(CreatePheGenIModulesFixture, CorrectProteinModulesMembers) {
    // Check correct members
    EXPECT_EQ(4, protein_modules.group_to_members[traits.str_to_int["GROUP1"]].count());
    // -- Check a gene with two protein products appears in the new modules
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP1"]][proteins.str_to_int["A1"]]);
    // -- Check that one of the protein products of a gene gets removed of the modules, because that one does not interact with other members
    EXPECT_FALSE(protein_modules.group_to_members[traits.str_to_int["GROUP1"]][proteins.str_to_int["A2"]]);
    // -- Check that one of the protein products of a gene gets removed of the modules, because that one does not interact with other members
    EXPECT_FALSE(protein_modules.group_to_members[traits.str_to_int["GROUP1"]][proteins.str_to_int["A3"]]);
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP1"]][proteins.str_to_int["B1"]]);
    // -- Check a gene with two protein products appears in the new modules
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP1"]][proteins.str_to_int["C1"]]);
    // -- Check that one of the protein products of a gene gets removed of the modules, because that one does not interact with other members
    EXPECT_FALSE(protein_modules.group_to_members[traits.str_to_int["GROUP1"]][proteins.str_to_int["C2"]]);
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP1"]][proteins.str_to_int["C3"]]);
    EXPECT_FALSE(protein_modules.group_to_members[traits.str_to_int["GROUP1"]][proteins.str_to_int["F1"]]);

    EXPECT_EQ(6, protein_modules.group_to_members[traits.str_to_int["GROUP2"]].count());
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP2"]][proteins.str_to_int["D1"]]);
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP2"]][proteins.str_to_int["D2"]]);
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP2"]][proteins.str_to_int["E1"]]);
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP2"]][proteins.str_to_int["F1"]]);
    EXPECT_FALSE(protein_modules.group_to_members[traits.str_to_int["GROUP2"]][proteins.str_to_int["F2"]]);
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP2"]][proteins.str_to_int["G1"]]);
    EXPECT_TRUE(protein_modules.group_to_members[traits.str_to_int["GROUP2"]][proteins.str_to_int["G2"]]);
}

TEST_F(CreatePheGenIModulesFixture, CorrectProteinModulesOwners) {
    // -- -- Check the proteins have/do not have ownership of the traits
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["A1"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["A2"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["A3"]].count());
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["B1"]].count());
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["C1"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["C2"]].count());
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["C3"]].count());

    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["A1"]][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["B1"]][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["C1"]][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["C3"]][traits.str_to_int["GROUP1"]]);

    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["D1"]].count());
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["D2"]].count());
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["E1"]].count());
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["F1"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["F2"]].count());
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["G1"]].count());
    EXPECT_EQ(1, protein_modules.member_to_groups[proteins.str_to_int["G2"]].count());

    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["D1"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["D2"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["E1"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["F1"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["G1"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(protein_modules.member_to_groups[proteins.str_to_int["G2"]][traits.str_to_int["GROUP2"]]);

    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["H1"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["H2"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["I1"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["J1"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["J2"]].count());
    EXPECT_EQ(0, protein_modules.member_to_groups[proteins.str_to_int["J3"]].count());
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
    EXPECT_EQ(2, proteoform_modules.group_to_members.size());
    EXPECT_EQ(2, proteoform_modules.member_to_groups[0].size());

    // Check there are 20 proteins
    EXPECT_EQ(26, proteoform_modules.member_to_groups.size());
    EXPECT_EQ(26, proteoform_modules.group_to_members[0].size());
}

TEST_F(CreatePheGenIModulesFixture, CorrectProteoformModulesMembers) {

    EXPECT_EQ(5, proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]].count());
    EXPECT_FALSE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["A1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["A1_2"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["A1_3"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["B1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["B1_2"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["C1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["C1_2"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["C2_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP1"]][proteoforms.str_to_int["C3_1"]]);

    EXPECT_EQ(8, proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]].count());
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["D1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["D1_2"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["D2_1"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["E1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["E1_2"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["F1_1"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["F2_1"]]);
    EXPECT_FALSE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["F2_2"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["G1_1"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["G1_2"]]);
    EXPECT_TRUE(proteoform_modules.group_to_members[traits.str_to_int["GROUP2"]][proteoforms.str_to_int["G2_1"]]);
}

TEST_F(CreatePheGenIModulesFixture, CorrectProtoformModulesOwners) {
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["A1_1"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["A1_2"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["A1_3"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["B1_1"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["B1_2"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["C1_1"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["C1_2"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["C2_1"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["C3_1"]].count());

    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["A1_2"]][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["B1_1"]][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["B1_2"]][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["C1_2"]][traits.str_to_int["GROUP1"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["C3_1"]][traits.str_to_int["GROUP1"]]);

    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["D1_1"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["D1_2"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["D2_1"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["E1_1"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["E1_2"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["F1_1"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["F2_1"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["F2_2"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["G1_1"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["G1_2"]].count());
    EXPECT_EQ(1, proteoform_modules.member_to_groups[proteoforms.str_to_int["G2_1"]].count());

    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["D1_1"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["D1_2"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["D2_1"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["E1_2"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["F1_1"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["G1_1"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["G1_2"]][traits.str_to_int["GROUP2"]]);
    EXPECT_TRUE(proteoform_modules.member_to_groups[proteoforms.str_to_int["G2_1"]][traits.str_to_int["GROUP2"]]);

    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["H1_1"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["H2_1"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["I1_1"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["J1_1"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["J2_1"]].count());
    EXPECT_EQ(0, proteoform_modules.member_to_groups[proteoforms.str_to_int["J3_1"]].count());
}

