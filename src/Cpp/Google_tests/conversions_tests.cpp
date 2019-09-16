#include <conversions.hpp>
#include "gtest/gtest.h"
#include "types.hpp"

class ConversionsFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        myStrSet.insert("A");
        myStrSet.insert("D");
        myStrSet.insert("B");
        myStrSet.insert("C");

        myStrVec = convert_uss_to_vs(myStrSet);
    }

    uss myStrSet;
    vs myStrVec;
};

TEST_F(ConversionsFixture, SameNumberOfElements) {
    EXPECT_EQ(myStrSet.size(), myStrVec.size());
}


TEST_F(ConversionsFixture, SameElements) {
    for (const auto &vec_element : myStrVec) {
        EXPECT_TRUE(myStrSet.find(vec_element) != myStrSet.end());
    }
}

class ReadMappingFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        mapping = readMapping(path_file_mapping);
        mapping_without_header = readMapping(path_file_mapping, false);
    }

    std::string path_file_mapping = "../../../../resources/UniProt/mapping_protein_to_genes_slice.tab";
    entity_mapping mapping;
    entity_mapping mapping_without_header;
};

// Test read mapping: Read file without header row
TEST_F(ReadMappingFixture, FileWithoutHeader){
    EXPECT_EQ(15, mapping.first_to_second.size()) << "There should be 16 genes";
    EXPECT_EQ(7, mapping.second_to_first.size()) << "There should be 7 proteins";
}

// Test read mapping: Correct number of sources, ignore header row by default
TEST_F(ReadMappingFixture, CorrectNumberOfSources) {
    EXPECT_EQ(13, mapping.first_to_second.size()) << "There should be 13 genes";
    EXPECT_EQ(6, mapping.second_to_first.size()) << "There should be 6 proteins";
}

// Test read mapping: Correct mapping for some sources
TEST_F(ReadMappingFixture, CorrectMapping) {
    std::set<std::string> values;
    auto ret = mapping.first_to_second.equal_range("P18825");
    for (auto it = ret.first; it != ret.second; ++it)
        values.insert(it->second);
    EXPECT_EQ(3, mapping.first_to_second.count("P18825"));
    EXPECT_TRUE(values.find("ADRA2C") != values.end());
    EXPECT_TRUE(values.find("ADRA2L2") != values.end());
    EXPECT_TRUE(values.find("ADRA2RL2") != values.end());
    values.clear();

    ret = mapping.first_to_second.equal_range("Q99424");
    for (auto it = ret.first; it != ret.second; ++it)
        values.insert(it->second);
    EXPECT_EQ(1, mapping.first_to_second.count("Q99424"));
    EXPECT_TRUE(values.find("ACOX2") != values.end());
    values.clear();

    ret = mapping.first_to_second.equal_range("Q8WW27");
    for (auto it = ret.first; it != ret.second; ++it)
        values.insert(it->second);
    EXPECT_EQ(2, mapping.first_to_second.count("Q8WW27"));
    EXPECT_TRUE(values.find("APOBEC4") != values.end());
    EXPECT_TRUE(values.find("C1orf169") != values.end());
    values.clear();

    EXPECT_FALSE(mapping.first_to_second.find("Entry") != mapping.first_to_second.end());
}

// Test read mapping: Correct total number of destinations
TEST_F(ReadMappingFixture, CorrectInverseMapping) {
    EXPECT_FALSE(mapping.second_to_first.find("Gene names") != mapping.first_to_second.end());

    std::set<std::string> values;
    auto ret = mapping.second_to_first.equal_range("ANKRD20A5P");
    for (auto it = ret.first; it != ret.second; ++it)
        values.insert(it->second);
    EXPECT_EQ(1, mapping.first_to_second.count("ANKRD20A5P"));
    EXPECT_TRUE(values.find("A0PJZ0") != values.end());
    values.clear();

    ret = mapping.second_to_first.equal_range("IDGFL");
    for (auto it = ret.first; it != ret.second; ++it)
        values.insert(it->second);
    EXPECT_EQ(1, mapping.first_to_second.count("IDGFL"));
    EXPECT_TRUE(values.find("Q9NZK5") != values.end());
    values.clear();
}
