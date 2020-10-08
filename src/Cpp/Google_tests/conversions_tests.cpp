#include "gtest/gtest.h"
#include <bimap_int_int.hpp>
#include "types.hpp"
#include "maps.hpp"

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

    std::string path_file_mapping = "../../Google_tests/resources/mapping_protein_to_genes_slice.tab";
    bimap_str_str mapping;
    bimap_str_str mapping_without_header;
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
    EXPECT_EQ(13, getValues(mapping.first_to_second).size()) << "There should be 15 genes";
    EXPECT_EQ(6, getKeys(mapping.first_to_second).size()) << "There should be 7 proteins";
    EXPECT_EQ(13, getKeys(mapping.second_to_first).size()) << "There should be 15 genes";
    EXPECT_EQ(6, getValues(mapping.second_to_first).size()) << "There should be 7 proteins";
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


class ReadMappingWithMultipleColumnsFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        mapping = readMapping(path_file_mapping, true, true);
    }

    bimap_str_str mapping;
    std::string path_file_mapping = "../../Google_tests/resources/createPheGenIModules/mapping_proteins_to_proteoforms.csv";
};


TEST_F(ReadMappingWithMultipleColumnsFixture, CorrectKeys) {

    EXPECT_EQ(18, getKeys(mapping.first_to_second).size());
    EXPECT_FALSE(hasKey(mapping.first_to_second, (std::string) "PROTEINS"));
    EXPECT_FALSE(hasKey(mapping.first_to_second, (std::string) "EXTRA_COLUMN"));
    EXPECT_FALSE(hasKey(mapping.first_to_second, (std::string) "extra"));
    EXPECT_TRUE(hasKey(mapping.first_to_second, (std::string) "A1"));
    EXPECT_TRUE(hasKey(mapping.first_to_second, (std::string) "C1"));
    EXPECT_TRUE(hasKey(mapping.first_to_second, (std::string) "H1"));
}

TEST_F(ReadMappingWithMultipleColumnsFixture, CorrectValues) {
    EXPECT_EQ(26, getValues(mapping.first_to_second).size());
    EXPECT_FALSE(hasValue(mapping.first_to_second, (std::string) "PROTEOFORMS"));
    EXPECT_FALSE(hasValue(mapping.first_to_second, (std::string) "EXTRA_COLUMN"));
    EXPECT_FALSE(hasValue(mapping.first_to_second, (std::string) "extra"));
    EXPECT_TRUE(hasValue(mapping.first_to_second, (std::string) "A1_1"));
    EXPECT_TRUE(hasValue(mapping.first_to_second, (std::string) "C1_2"));
    EXPECT_TRUE(hasValue(mapping.first_to_second, (std::string) "J3_1"));
}

TEST_F(ReadMappingWithMultipleColumnsFixture, CorrectInverseKeys) {
    EXPECT_EQ(26, getKeys(mapping.second_to_first).size());
    EXPECT_FALSE(hasKey(mapping.second_to_first, (std::string) "PROTEOFORMS"));
    EXPECT_FALSE(hasKey(mapping.second_to_first, (std::string) "EXTRA_COLUMN"));
    EXPECT_FALSE(hasKey(mapping.second_to_first, (std::string) "extra"));
    EXPECT_TRUE(hasKey(mapping.second_to_first, (std::string) "A1_1"));
    EXPECT_TRUE(hasKey(mapping.second_to_first, (std::string) "C1_2"));
    EXPECT_TRUE(hasKey(mapping.second_to_first, (std::string) "J3_1"));
}

TEST_F(ReadMappingWithMultipleColumnsFixture, CorrectInverseValues) {
    EXPECT_EQ(18, getValues(mapping.second_to_first).size());
    EXPECT_FALSE(hasValue(mapping.second_to_first, (std::string) "PROTEINS"));
    EXPECT_FALSE(hasValue(mapping.second_to_first, (std::string) "EXTRA_COLUMN"));
    EXPECT_FALSE(hasValue(mapping.second_to_first, (std::string) "extra"));
    EXPECT_TRUE(hasValue(mapping.second_to_first, (std::string) "A1"));
    EXPECT_TRUE(hasValue(mapping.second_to_first, (std::string) "C1"));
    EXPECT_TRUE(hasValue(mapping.second_to_first, (std::string) "H1"));
}


class ReadMappingWithMultipleColumnsFileHasNoHeaderFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        mapping = readMapping(path_file_mapping, false, true);
    }

    bimap_str_str mapping;
    std::string path_file_mapping = "../../Google_tests/resources/createPheGenIModules/mapping_proteins_to_proteoforms.csv";
};

TEST_F(ReadMappingWithMultipleColumnsFileHasNoHeaderFixture, FileWithoutHeaderCorrectKeys) {

    EXPECT_EQ(19, getKeys(mapping.first_to_second).size());
    EXPECT_TRUE(hasKey(mapping.first_to_second, (std::string) "PROTEINS"));
    EXPECT_FALSE(hasKey(mapping.first_to_second, (std::string) "EXTRA_COLUMN"));
    EXPECT_FALSE(hasKey(mapping.first_to_second, (std::string) "extra"));
    EXPECT_TRUE(hasKey(mapping.first_to_second, (std::string) "A1"));
    EXPECT_TRUE(hasKey(mapping.first_to_second, (std::string) "J3"));
}

TEST_F(ReadMappingWithMultipleColumnsFileHasNoHeaderFixture, FileWithoutHeaderCorrectValues) {
    EXPECT_EQ(27, getValues(mapping.first_to_second).size());
    EXPECT_TRUE(hasValue(mapping.first_to_second, (std::string) "PROTEOFORMS"));
    EXPECT_FALSE(hasValue(mapping.first_to_second, (std::string) "EXTRA_COLUMN"));
    EXPECT_FALSE(hasValue(mapping.first_to_second, (std::string) "extra"));
    EXPECT_TRUE(hasValue(mapping.first_to_second, (std::string) "A1_1"));
    EXPECT_TRUE(hasValue(mapping.first_to_second, (std::string) "C1_2"));
    EXPECT_TRUE(hasValue(mapping.first_to_second, (std::string) "J3_1"));
}

TEST_F(ReadMappingWithMultipleColumnsFileHasNoHeaderFixture, CorrectInverseKeys) {
    EXPECT_EQ(27, getKeys(mapping.second_to_first).size());
    EXPECT_TRUE(hasKey(mapping.second_to_first, (std::string) "PROTEOFORMS"));
    EXPECT_FALSE(hasKey(mapping.second_to_first, (std::string) "EXTRA_COLUMN"));
    EXPECT_FALSE(hasKey(mapping.second_to_first, (std::string) "extra"));
    EXPECT_TRUE(hasKey(mapping.second_to_first, (std::string) "A1_1"));
    EXPECT_TRUE(hasKey(mapping.second_to_first, (std::string) "C1_2"));
    EXPECT_TRUE(hasKey(mapping.second_to_first, (std::string) "J3_1"));
}

TEST_F(ReadMappingWithMultipleColumnsFileHasNoHeaderFixture, CorrectInverseValues) {
    EXPECT_EQ(19, getValues(mapping.second_to_first).size());
    EXPECT_TRUE(hasValue(mapping.second_to_first, (std::string) "PROTEINS"));
    EXPECT_FALSE(hasValue(mapping.second_to_first, (std::string) "EXTRA_COLUMN"));
    EXPECT_FALSE(hasValue(mapping.second_to_first, (std::string) "extra"));
    EXPECT_TRUE(hasValue(mapping.second_to_first, (std::string) "A1"));
    EXPECT_TRUE(hasValue(mapping.second_to_first, (std::string) "C1"));
    EXPECT_TRUE(hasValue(mapping.second_to_first, (std::string) "H1"));
}