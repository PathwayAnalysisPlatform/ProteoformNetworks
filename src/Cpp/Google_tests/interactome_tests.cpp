#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <vector>
#include <Interactome.hpp>
#include <tuple>

using ::testing::UnorderedElementsAreArray;
using ::testing::Contains;
using ::testing::AnyOf;
using ::testing::Not;

class InteractomeFixture : public ::testing::Test {
private:
    std::vector<std::pair<int, int>> interactions = {
            std::make_pair(1, 2),
            std::make_pair(2, 3),
            std::make_pair(4, 5)
    };
protected:
    Interactome interactome;

    InteractomeFixture() : interactome{interactions} {}

    virtual void SetUp() override {
        std::string names = "1 A\n2 B\n3 C\n4 D\n5 E";
        std::istringstream ss(names);
        interactome.readNodeNames(ss);
    }

    virtual void TearDown() override {

    }
};

TEST_F(InteractomeFixture, CanCallConstructorWithRangesAndInteractionsAndNames){
    // TODO
}

TEST_F(InteractomeFixture, GetNodesAsVectorTest) {
    ASSERT_EQ(typeid(interactome.getNodes()), typeid(std::vector<int>));
    ASSERT_THAT(interactome.getNodes(), UnorderedElementsAreArray({1, 2, 3, 4, 5}));
}

TEST_F(InteractomeFixture, AddNodeTest) {
    interactome.addNode(7);
    ASSERT_THAT(interactome.getNodes(), Contains(7));
}

TEST_F(InteractomeFixture, CreateNetworkFromInteractionList) {
    std::vector<std::pair<int, int>> interactions = {
            std::make_pair(11, 12),
            std::make_pair(12, 13),
            std::make_pair(14, 15)
    };
    interactome = Interactome(interactions);

    ASSERT_EQ(5, interactome.getNodes().size());
    ASSERT_THAT(interactome.getNodes(), UnorderedElementsAreArray({11, 12, 13, 14, 15}));
}

TEST_F(InteractomeFixture, GetAllInteractorsTest) {
    std::vector<int> neighbors = interactome.getInteractors(2);

    ASSERT_EQ(neighbors.size(), 2);
    ASSERT_THAT(neighbors, UnorderedElementsAreArray({1, 3}));
    ASSERT_THAT(neighbors, Not(AnyOf(Contains(4), Contains(5))));
}

TEST_F(InteractomeFixture, GetAllInteractorsWithNonexistentNodeTest) {
    std::vector<int> neighbors;
    ASSERT_THROW(neighbors = interactome.getInteractors(7), std::out_of_range);
}

TEST_F(InteractomeFixture, CheckIfNodeExistsByIndexTest) {
    ASSERT_TRUE(interactome.hasNode(1));
    ASSERT_FALSE(interactome.hasNode(7));
}

TEST_F(InteractomeFixture, GetNodeNamesTest) {
    std::string name;
    ASSERT_EQ(name = interactome.getNodeName(1), "A");
    ASSERT_EQ(name = interactome.getNodeName(5), "E");
    ASSERT_THROW(name = interactome.getNodeName(7), std::invalid_argument);
}

TEST_F(InteractomeFixture, ReadNodeNamesTest) {
    // The names were read at SetUp
    std::string names = "1 A2\n3 C2\n2 B2\n4 D2\n5 E2";
    std::istringstream ss(names);
    interactome.readNodeNames(ss);
    ASSERT_EQ(interactome.getNodeName(2), "B2");
    ASSERT_EQ(interactome.getNodeName(5), "E2");
}

TEST_F(InteractomeFixture, ReadNodeNamesAgainReplacesNamesTest) {
    std::string names = "1 A1\n3 C1\n2 B1\n4 D1\n5 E1";
    std::istringstream ss(names);
    interactome.readNodeNames(ss);
    names = "1 A2\n3 C2\n2 B2\n4 D2\n5 E2";
    ss = std::istringstream(names);
    interactome.readNodeNames(ss);
    ASSERT_EQ(interactome.getNodeName(2), "B2");
    ASSERT_EQ(interactome.getNodeName(5), "E2");
}

TEST_F(InteractomeFixture, RepeatedNodeNamesThrowsExceptionTest) {
    std::string names = "1 A\n2 B\n3 A\n4 D\n5 E";
    std::istringstream ss(names);
    ASSERT_THROW(interactome.readNodeNames(ss), std::invalid_argument);
}

TEST_F(InteractomeFixture, RepeatedNodesInTheNamesStreamThrowsException) {
    std::string names = "1 A\n2 B\n2 C\n4 D\n5 E";
    std::istringstream ss(names);
    ASSERT_THROW(interactome.readNodeNames(ss), std::invalid_argument);
}

TEST_F(InteractomeFixture, InitialNodeNamesAreSetTest) {
    std::vector<std::pair<int, int>> interactions1({std::make_pair(1, 2)});
    Interactome interactome1(interactions1);
    std::string name;
    ASSERT_EQ(name = interactome1.getNodeName(1), "1");
    ASSERT_EQ(name = interactome1.getNodeName(2), "2");
}

TEST_F(InteractomeFixture, ThrowExceptionWhenRequestNameNonExistentNodeTest) {
    std::string name;
    ASSERT_THROW(name = interactome.getNodeName(8), std::invalid_argument);
}

TEST_F(InteractomeFixture, NamesStreamHasLessNamesThrowsExceptionTest) {
    std::string names = "1 A\n2 B\n3 C\n4 D";
    std::istringstream ss(names);
    ASSERT_THROW(interactome.readNodeNames(ss), std::invalid_argument);
}

TEST_F(InteractomeFixture, ThrowsExceptionWithNodeNotInNodesTest) {
    std::string names = "6 F\n1 A";
    std::stringstream ss(names);
    try {
        interactome.readNodeNames(ss);
        FAIL();
    } catch (std::invalid_argument const &e) {
        ASSERT_EQ(e.what(), std::string("Provided name for unexistent node: 6"));
    }
}

TEST_F(InteractomeFixture, NamesStreamHasMoreNamesThanExpectedTest) {
    std::string names = "1 A\n2 B\n3 C\n4 D\n5 E\n6 F\n";
    std::istringstream ss(names);
    try {
        interactome.readNodeNames(ss);
        FAIL();
    } catch (std::invalid_argument &e) {
        ASSERT_EQ(e.what(), std::string("Provided too many arguments to name the nodes."));
    }
}

TEST_F(InteractomeFixture, CheckIfNodeExistsByNameTest) {
    std::string name = "A";
    ASSERT_TRUE(interactome.hasNode(name));
    ASSERT_FALSE(interactome.hasNode("F"));
}

TEST_F(InteractomeFixture, ReadTypesTest){
    std::string ranges = "";
    std::istringstream ss(ranges);
    interactome.readTypeRanges(ss);

    ASSERT_EQ(interactome.getGenesStart(0));
}