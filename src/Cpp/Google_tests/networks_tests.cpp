#include <bimap_str_int.hpp>
#include <networks.hpp>
#include "gtest/gtest.h"
#include "types.hpp"
#include "maps.hpp"

class LoadInteractionNetworkFixture : public testing::Test {
protected:
    virtual void SetUp() {
        proteins = createBimap(path_file_proteins);
        interactions = loadInteractionNetwork(path_file_interactions, proteins, true);
    }

    ummii interactions;
    std::string path_file_proteins = "../../Google_tests/resources/proteins_slice.csv";
    bimap_str_int proteins;
    std::string path_file_interactions = "../../Google_tests/resources/proteinInternalEdges_slice.tsv";
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
//        std::cerr << "The execution path is: " << get_current_dir_name() << std::endl;

        auto ret = loadModules(path_file_modules);
        example_modules = ret.entity_modules;
        groups = ret.groups;
        members = ret.members;
    }

    std::string path_file_modules = "../../Google_tests/resources/example_modules.csv";
    modules example_modules;
    bimap_str_int groups, members;
};

// Check Members correct size
TEST_F(LoadModulesFixture, CorrectMembersSize) {
    EXPECT_EQ(5, members.str_to_int.size());
    EXPECT_EQ(5, members.int_to_str.size());
}

// Check Members int_to_str elements are correct
TEST_F(LoadModulesFixture, CorrectMembers) {
    EXPECT_EQ("1", members.int_to_str[0]);
    EXPECT_EQ("3", members.int_to_str[2]);
    EXPECT_EQ("5", members.int_to_str[4]);
}

// Check Members str_to_int keys are correct
TEST_F(LoadModulesFixture, CorrectMembersStrToIntKeys) {
    ASSERT_TRUE(hasKey(members.str_to_int, "1"));
    ASSERT_TRUE(hasKey(members.str_to_int, "3"));
    ASSERT_TRUE(hasKey(members.str_to_int, "5"));
}

// Check Members str_to_int index values are correct
TEST_F(LoadModulesFixture, CorrectMembersStrToIntValues) {
    EXPECT_EQ(0, members.str_to_int["1"]);
    EXPECT_EQ(2, members.str_to_int["3"]);
    EXPECT_EQ(4, members.str_to_int["5"]);
}

// Check groups correct size
TEST_F(LoadModulesFixture, CorrectGroupsSize) {
    EXPECT_EQ(3, groups.str_to_int.size());
    EXPECT_EQ(3, groups.int_to_str.size());
}

// Check Groups int_to_str elements are correct
TEST_F(LoadModulesFixture, CorrectGroups) {
    EXPECT_EQ("A", groups.int_to_str[0]);
    EXPECT_EQ("B", groups.int_to_str[1]);
    EXPECT_EQ("C", groups.int_to_str[2]);
}

// Check Groups str_to_int keys are correct
TEST_F(LoadModulesFixture, CorrectGroupsStrToIntKeys) {
    ASSERT_TRUE(hasKey(groups.str_to_int, "A"));
    ASSERT_TRUE(hasKey(groups.str_to_int, "B"));
    ASSERT_TRUE(hasKey(groups.str_to_int, "C"));
}

// Check Groups str_to_int index values are correct
TEST_F(LoadModulesFixture, CorrectGroupsStrToIntValues) {
    EXPECT_EQ(0, groups.str_to_int["A"]);
    EXPECT_EQ(1, groups.str_to_int["B"]);
    EXPECT_EQ(2, groups.str_to_int["C"]);
}

// Check modules sizes
TEST_F(LoadModulesFixture, CorrectModulesSizes) {
    EXPECT_EQ(3, example_modules.group_to_members.size());
    for (auto group_entry = example_modules.group_to_members.begin();
         group_entry != example_modules.group_to_members.end(); group_entry++) {
        EXPECT_EQ(groups.int_to_str.size(), group_entry->second.size());
    }
    EXPECT_EQ(5, example_modules.member_to_groups.size());
    for (auto member_entry = example_modules.member_to_groups.begin();
         member_entry != example_modules.member_to_groups.end(); member_entry++) {
        EXPECT_EQ(members.int_to_str.size(), member_entry->second.size());
    }
}

// Check modules group --> members
TEST_F(LoadModulesFixture, CorrectGroupToMembers) {
    // Correct number of members
    EXPECT_EQ(3, example_modules.group_to_members["A"].count());
    EXPECT_EQ(2, example_modules.group_to_members["B"].count());
    EXPECT_EQ(1, example_modules.group_to_members["C"].count());
    // Correct members
    EXPECT_TRUE(example_modules.group_to_members["A"][members.str_to_int["1"]]);
    EXPECT_TRUE(example_modules.group_to_members["A"][members.str_to_int["2"]]);
    EXPECT_TRUE(example_modules.group_to_members["A"][members.str_to_int["3"]]);

    EXPECT_TRUE(example_modules.group_to_members["B"][members.str_to_int["3"]]);
    EXPECT_TRUE(example_modules.group_to_members["B"][members.str_to_int["4"]]);

    EXPECT_TRUE(example_modules.group_to_members["C"][members.str_to_int["5"]]);
}

// Check modules member --> groups
TEST_F(LoadModulesFixture, CorrectMemberToGroups) {
    // Correct number of owners
    EXPECT_EQ(1, example_modules.member_to_groups["1"].count());
    EXPECT_EQ(1, example_modules.member_to_groups["2"].count());
    EXPECT_EQ(2, example_modules.member_to_groups["3"].count());
    EXPECT_EQ(1, example_modules.member_to_groups["4"].count());
    EXPECT_EQ(1, example_modules.member_to_groups["5"].count());

    // Correct owners
    EXPECT_TRUE(example_modules.member_to_groups["1"][groups.str_to_int["A"]]);
    EXPECT_TRUE(example_modules.member_to_groups["2"][groups.str_to_int["A"]]);
    EXPECT_TRUE(example_modules.member_to_groups["3"][groups.str_to_int["A"]]);
    EXPECT_TRUE(example_modules.member_to_groups["3"][groups.str_to_int["B"]]);
    EXPECT_TRUE(example_modules.member_to_groups["4"][groups.str_to_int["B"]]);
    EXPECT_TRUE(example_modules.member_to_groups["5"][groups.str_to_int["C"]]);
}

class RemoveDisconnectedMembersFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        auto[example_modules, groups, members] = loadModules(path_file_modules);
        interactions = loadInteractionNetwork(path_file_interactions, members, true);
    }

    std::string path_file_modules = "../../Google_tests/resources/example_modules.csv";
    std::string path_file_interactions = "../../Google_tests/resources/example_interactions.csv";
    ummii interactions;

};


TEST_F(RemoveDisconnectedMembersFixture, KeepsConnectedMembers) {
    auto[example_modules, groups, members] = loadModules(path_file_modules);

    // Check vertices are members
    ASSERT_TRUE(example_modules.group_to_members["A"][members.str_to_int["1"]]);
    ASSERT_TRUE(example_modules.group_to_members["A"][members.str_to_int["2"]]);
    ASSERT_TRUE(example_modules.group_to_members["A"][members.str_to_int["3"]]);

    ASSERT_TRUE(example_modules.group_to_members["B"][members.str_to_int["3"]]);
    ASSERT_TRUE(example_modules.group_to_members["B"][members.str_to_int["4"]]);

    removeDisconnectedMembers(example_modules, groups, members, interactions);

    // In module A, should keep 1 and 2. But delete vertex 3, because it is disconnected
    EXPECT_TRUE(example_modules.group_to_members["A"][members.str_to_int["1"]]);
    EXPECT_TRUE(example_modules.group_to_members["A"][members.str_to_int["2"]]);

    // In module B, should keep 3 and 4
    EXPECT_TRUE(example_modules.group_to_members["B"][members.str_to_int["3"]]);
    EXPECT_TRUE(example_modules.group_to_members["B"][members.str_to_int["4"]]);
}

TEST_F(RemoveDisconnectedMembersFixture, RemovesDisconnectedMember) {
    auto[example_modules, groups, members] = loadModules(path_file_modules);

    // Check vertices are members
    ASSERT_TRUE(example_modules.group_to_members["A"][members.str_to_int["3"]]);

    ASSERT_TRUE(example_modules.group_to_members["C"][members.str_to_int["5"]]);

    auto ret = removeDisconnectedMembers(example_modules, groups, members, interactions);

    // In module A, should keep 1 and 2. But delete vertex 3, because it is disconnected
    // Removes vertex 3 from module A
    EXPECT_FALSE(example_modules.group_to_members["A"][members.str_to_int["3"]]);

    // Removes vertex 5 from module C
    // In module C, should remove the only member 5
    EXPECT_FALSE(example_modules.group_to_members["C"][members.str_to_int["5"]]);
}
