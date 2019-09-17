#include "gtest/gtest.h"
#include <conversions.hpp>
#include "types.hpp"

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
TEST_F(ReadMappingFixture, FileWithoutHeader) {
    std::set<std::string> proteins, genes;
    for (auto it = mapping_without_header.first_to_second.begin();
         it != mapping_without_header.first_to_second.end();
         it++) {
        proteins.insert(it->first);
        genes.insert(it->second);
    }
    EXPECT_EQ(15, genes.size()) << "There should be 15 genes";
    EXPECT_EQ(7, proteins.size()) << "There should be 7 proteins";
    EXPECT_EQ(15, mapping_without_header.second_to_first.size()) << "There should be 15 map entries";
    EXPECT_EQ(15, mapping_without_header.second_to_first.size()) << "There should be 15 map entries";
}

// Test read mapping: Correct number of sources, ignore header row by default
TEST_F(ReadMappingFixture, CorrectNumberOfSources) {
    std::set<std::string> proteins, genes;
    for (auto it = mapping.first_to_second.begin(); it != mapping.first_to_second.end(); it++) {
        proteins.insert(it->first);
        genes.insert(it->second);
    }
    EXPECT_EQ(13, genes.size()) << "There should be 15 genes";
    EXPECT_EQ(6, proteins.size()) << "There should be 7 proteins";
    EXPECT_EQ(13, mapping.first_to_second.size()) << "There should be 13 map entries";
    EXPECT_EQ(13, mapping.second_to_first.size()) << "There should be 13 map entries";
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
    EXPECT_EQ(1, mapping.first_to_second.count("Q8WW27"));
    EXPECT_TRUE(values.find("APOBEC4") != values.end());
    values.clear();

    EXPECT_FALSE(mapping.first_to_second.find("Entry") != mapping.first_to_second.end());
}

// Test read mapping: Correct total number of destinations
TEST_F(ReadMappingFixture, CorrectInverseMapping) {
    EXPECT_FALSE(mapping.second_to_first.find("Gene names") != mapping.second_to_first.end());

    EXPECT_EQ(1, mapping.second_to_first.count("ANKRD20A5P"));
    EXPECT_EQ("A0PJZ0", (mapping.second_to_first.equal_range("ANKRD20A5P")).first->second);

    std::set<std::string> values;
    auto ret = mapping.second_to_first.equal_range("IDGFL");
    for (auto it = ret.first; it != ret.second; ++it)
        values.insert(it->second);
    EXPECT_EQ(1, mapping.second_to_first.count("IDGFL"));
    EXPECT_TRUE(values.find("Q9NZK5") != values.end());
    values.clear();
}