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
    std::string path_file_proteins = "../../../../resources/Reactome/v70/Proteins/proteins_slice.csv";
    bimap_str_int proteins;
    std::string path_file_interactions = "../../../../resources/Reactome/v70/Proteins/proteinInternalEdges_slice.tsv";
};


TEST_F(LoadInteractionNetworkFixture, CorrectTotalEntities) {
    EXPECT_EQ(300, interactions.size()) << "Expected double number of interactions in the file, due to redundant inverse edges.";
}


TEST_F(LoadInteractionNetworkFixture, CorrectNeighbours) {
    ASSERT_TRUE(interactions.find(proteins.str_to_int["P19224"]) != interactions.end());
    auto range = interactions.equal_range(proteins.str_to_int["P19224"]);
    std::set<std::string> neighbors;
    for(auto it = range.first; it != range.second; it++){
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

TEST_F(LoadInteractionNetworkFixture, InverseEdges){
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
// TODO: Test reading interactions from file