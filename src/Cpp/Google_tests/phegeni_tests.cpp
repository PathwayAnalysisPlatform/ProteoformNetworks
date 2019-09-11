#include "gtest/gtest.h"
#include <bimap_str_int.hpp>
#include "phegeni.hpp"

#include <windows.h>

std::string GetExeFileName()
{
    char buffer[MAX_PATH];
    GetModuleFileName( NULL, buffer, MAX_PATH );
    return std::string(buffer);
}

std::string GetExePath()
{
    std::string f = GetExeFileName();
    return f.substr(0, f.find_last_of( "\\/" ));
}

class PhegeniLoadGenesAndTraitsFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        std::cout << GetExePath() << std::endl;
        path_file_phegeni = "../../../../resources/PheGenI/PheGenI_Association_genome_wide_significant_slice.txt";
        path_file_genes = "../../../../resources/Reactome/v69/genes_slice.csv";
    }

    std::string path_file_phegeni, path_file_genes;

};
// Throws runtime error if file PheGenI is not found


TEST_F(PhegeniLoadGenesAndTraitsFixture, UnexistentFileRaisesError) {
    try {
        auto genesAndTraits = loadPheGenIGenesAndTraits("wrong_path", createBimap(path_file_genes));
        FAIL() << "Expected std::runtime_error";
    }
    catch(std::runtime_error const & err) {
        EXPECT_EQ(err.what(),std::string("Cannot open path_file_phegeni at loadPheGenIGenesAndTraits"));
    }
    catch(...) {
        FAIL() << "Expected std::runtime_error";
    }
}
// Check that genes not in Gene list are not added

TEST_F(PhegeniLoadGenesAndTraitsFixture, LoadPheGenIGenesAndTraits) {
    const bimap_str_int &genes = createBimap(path_file_genes);
    auto genesAndTraits = loadPheGenIGenesAndTraits(path_file_phegeni, genes);

    // Check that the traits are correct
    EXPECT_EQ(8, genesAndTraits.phegeni_traits.str_to_int.size());
    EXPECT_EQ(8, genesAndTraits.phegeni_traits.int_to_str.size());

    // Check that the genes are correct
    EXPECT_EQ(25, genesAndTraits.phegeni_genes.int_to_str.size());
    EXPECT_EQ(25, genesAndTraits.phegeni_genes.int_to_str.size());

    // Check genes are sorted
    EXPECT_EQ(1, genesAndTraits.phegeni_traits.str_to_int.at("\"Cholesterol, LDL\""));
    EXPECT_EQ("\"Cholesterol, LDL\"", genesAndTraits.phegeni_traits.int_to_str[1]);
    EXPECT_EQ(2, genesAndTraits.phegeni_traits.str_to_int.at("Bilirubin"));
    EXPECT_EQ("Bilirubin", genesAndTraits.phegeni_traits.int_to_str[2]);
    EXPECT_EQ(7, genesAndTraits.phegeni_traits.str_to_int.at("Wet Macular Degeneration"));
    EXPECT_EQ("Wet Macular Degeneration", genesAndTraits.phegeni_traits.int_to_str[7]);
}
