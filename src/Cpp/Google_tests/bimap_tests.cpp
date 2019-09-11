#include "gtest/gtest.h"
#include <iostream>
#include <list>
#include <bimap_str_int.hpp>

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

class BimapStrIntFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        std::cout << GetExePath() << std::endl;

        genes.push_back("MAP1ALC3");
        genes.push_back("PARK2");
        genes.push_back("PRKN");
        genes.push_back("PINK1");
        genes.push_back("TOMM40");
        genes.push_back("C19orf1");

        // Create artificial file
        std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
        file_name = test_name + "_list_file.csv";
        std::ofstream list_file(file_name);
        for (const auto &gene : genes)
            list_file << gene << "\n";
        list_file.close();
    }

    virtual void TearDown() {
        // Remove artificial file
        int n = file_name.length();
        char file_name_char[n + 1];
        strcpy(file_name_char, file_name.c_str());
        std::remove(file_name_char);
    }

    std::list<std::string> genes;
    std::string file_name;
};

TEST_F(BimapStrIntFixture, CreateBimap) {
    bimap_str_int bimap = createBimap(file_name, true);

    // Check right number of elements
    EXPECT_EQ(bimap.int_to_str.size(), 5);
    EXPECT_EQ(bimap.str_to_int.size(), 5);

    // Check the correct sorted elements
    EXPECT_EQ("C19orf1", bimap.int_to_str[0]);
    EXPECT_EQ("TOMM40", bimap.int_to_str[4]);

    EXPECT_EQ(0, bimap.str_to_int["C19orf1"]);
    EXPECT_EQ(4, bimap.str_to_int["TOMM40"]);
}

TEST_F(BimapStrIntFixture, CreateBimapWithoutHeader) {
    bimap_str_int bimap = createBimap(file_name, false);

    // Check right number of elements
    EXPECT_EQ(bimap.int_to_str.size(), 6);
    EXPECT_EQ(bimap.str_to_int.size(), 6);

    EXPECT_EQ("C19orf1", bimap.int_to_str[0]);
    EXPECT_EQ("TOMM40", bimap.int_to_str[5]);

    EXPECT_EQ(0, bimap.str_to_int["C19orf1"]);
    EXPECT_EQ(5, bimap.str_to_int["TOMM40"]);
}
