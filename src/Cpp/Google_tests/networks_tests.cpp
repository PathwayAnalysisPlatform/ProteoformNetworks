#include <bimap_str_int.hpp>
#include <networks.hpp>
#include "gtest/gtest.h"
#include "types.hpp"
#include "maps.hpp"

#include <windows.h>

std::string GetExeFileName() {
    char buffer[MAX_PATH];
    GetModuleFileName(NULL, buffer, MAX_PATH);
    return std::string(buffer);
}

std::string GetExePath() {
    std::string f = GetExeFileName();
    return f.substr(0, f.find_last_of("\\/"));
}

class LoadInteractionNetworkFixture : public testing::Test {
protected:
    virtual void SetUp() {
        proteins = createBimap(path_file_proteins);
        interactions = loadInteractionNetwork(path_file_interactions, proteins, true);
    }

    ummii interactions;
    std::string path_file_proteins = "../../../../resources/Reactome/v70/Proteins/proteins_slice.csv";
    bimap_str_int proteins;
    std::string path_file_interactions = "../../../../resources/Reactome/v70/Proteins/proteinInternalEdges_slice.tsv";
};


TEST_F(LoadInteractionNetworkFixture, CorrectTotalEntities) {
    EXPECT_EQ(300, interactions.size())
                        << "Expected double number of interactions in the file, due to redundant inverse edges.";
}


TEST_F(LoadInteractionNetworkFixture, CorrectNeighbours) {
    ASSERT_TRUE(interactions.find(proteins.str_to_int["P19224"]) != interactions.end());
    auto range = interactions.equal_range(proteins.str_to_int["P19224"]);
    std::set<std::string> neighbors;
    for (auto it = range.first; it != range.second; it++) {
        neighbors.insert(proteins.int_to_str[it->second]);
    }
    EXPECT_EQ(7, neighbors.size());
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P19224"], proteins.str_to_int["Q9HAW9"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P19224"], proteins.str_to_int["Q9HAW7"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P19224"], proteins.str_to_int["P35504"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P19224"], proteins.str_to_int["P35503"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P19224"], proteins.str_to_int["P22310"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P19224"], proteins.str_to_int["P22309"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P19224"], proteins.str_to_int["O60656"]));
}

TEST_F(LoadInteractionNetworkFixture, InverseEdges) {
    ASSERT_TRUE(interactions.find(proteins.str_to_int["Q9HAW9"]) != interactions.end());
    ASSERT_TRUE(interactions.find(proteins.str_to_int["Q9HAW7"]) != interactions.end());
    ASSERT_TRUE(interactions.find(proteins.str_to_int["P35504"]) != interactions.end());
    ASSERT_TRUE(interactions.find(proteins.str_to_int["P35503"]) != interactions.end());
    ASSERT_TRUE(interactions.find(proteins.str_to_int["P22310"]) != interactions.end());
    ASSERT_TRUE(interactions.find(proteins.str_to_int["P22309"]) != interactions.end());
    ASSERT_TRUE(interactions.find(proteins.str_to_int["O60656"]) != interactions.end());

    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["Q9HAW9"], proteins.str_to_int["P19224"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["Q9HAW7"], proteins.str_to_int["P19224"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P35504"], proteins.str_to_int["P19224"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P35503"], proteins.str_to_int["P19224"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P22310"], proteins.str_to_int["P19224"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["P22309"], proteins.str_to_int["P19224"]));
    EXPECT_TRUE(keyHasValue(interactions, proteins.str_to_int["O60656"], proteins.str_to_int["P19224"]));
}


class LoadModulesFixture : public ::testing::Test {
protected:
    virtual void SetUp() {
        std::cerr << GetExePath() << std::endl;
        auto ret = loadModules(path_file_modules);
        modules = ret.entity_modules;
        groups = ret.groups;
        members = ret.members;
    }

    std::string path_file_modules = "../../resources/modules.csv";
    modules modules;
    bimap_str_int groups, members;
};

TEST_F(LoadModulesFixture, CorrectMembersSize) {
    EXPECT_EQ(5, members.str_to_int.size());
    EXPECT_EQ(5, members.int_to_str.size());
}

// Check int_to_str elements are correct
TEST_F(LoadModulesFixture, CorrectMembers) {
    EXPECT_EQ("1", members.int_to_str[0]);
    EXPECT_EQ("3", members.int_to_str[2]);
    EXPECT_EQ("5", members.int_to_str[4]);
}

// Check str_to_int keys are correct
TEST_F(LoadModulesFixture, CorrectMembersStrToIntKeys){
    ASSERT_TRUE(hasKey(members.str_to_int, "A"));
    ASSERT_TRUE(hasKey(members.str_to_int, "B"));
    ASSERT_TRUE(hasKey(members.str_to_int, "C"));
}

// Check str_to_int values are correct
TEST_F(LoadModulesFixture, CorrectMembersStrToIntValues){
    EXPECT_EQ(0, members.str_to_int["A"]);
    EXPECT_EQ(1, members.str_to_int["B"]);
    EXPECT_EQ(2, members.str_to_int["C"]);
}

TEST_F(LoadModulesFixture, CorrectGroups) {
    // Check int_to_str size is correct

    // Check int_to_str elements are correct

    // Check str_to_int size is correct

    // Check str_to_int keys are correct

    // Check str_to_int values are correct
}

TEST_F(LoadModulesFixture, CorrectModules) {
    // Check keys in group_to_members are correct

    // Check values in group_to_members are correct

    // Check keys in member_to_groups are correct

    // Check values in member_to_groups are correct
}

class RemoveDisconnectedMembersFixture : public ::testing::Test {

protected:
    virtual void SetUp() {

        ummii interactions;
        removeDisconnectedMembers(example_modules, interactions);
    }

    modules example_modules;
};


TEST_F(RemoveDisconnectedMembersFixture, KeepsConnectedMembers) {

}

TEST_F(RemoveDisconnectedMembersFixture, RemovesDisconnectedMember) {
    // Check removes member from the members set
    // Check removes the set from the owners set
}
