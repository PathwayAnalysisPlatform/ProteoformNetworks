#include <bitset>
#include "gtest/gtest.h"
#include "scores.hpp"

class ScoresFixture : public ::testing::Test {

};

// Jaccard similarity when both sets are empty
TEST(ScoresFixture, JaccardBothEmpty) {
    base::dynamic_bitset<> set1, set2;

    EXPECT_EQ(0, set1.count()) << "The set1 should be empty.";
    EXPECT_EQ(0, set2.count()) << "The set2 should be empty.";
    EXPECT_EQ(1, getJaccardSimilarity(set1, set2));
    EXPECT_EQ(1, getJaccardSimilarity(set2, set1));

}

// Jaccard similarity when the sets are disjoint (completely different)
TEST(ScoresFixture, JaccardDisjoint) {
    base::dynamic_bitset<> set1(5), set2(5);
    set1[0].set();
    set1[1].set();
    set2[2].set();
    set2[3].set();
    set2[4].set();

    EXPECT_EQ(2, set1.count()) << "Set1 should have 2 members.";
    EXPECT_EQ(3, set2.count()) << "Set2 should have 3 members.";
    EXPECT_EQ(0, getJaccardSimilarity(set1, set2));
    EXPECT_EQ(0, getJaccardSimilarity(set2, set1));
}

// Jaccard similarity when both sets are the same
TEST(ScoresFixture, JaccardBothSame) {
    base::dynamic_bitset<> set1(5), set2(5);

    set1[0].set();
    set1[2].set();
    set1[4].set();
    set2[0].set();
    set2[2].set();
    set2[4].set();

    EXPECT_EQ(3, set1.count()) << "Set1 should have 3 members.";
    EXPECT_EQ(3, set2.count()) << "Set2 should have 3 members.";
    EXPECT_EQ(1, getJaccardSimilarity(set1, set2));
    EXPECT_EQ(1, getJaccardSimilarity(set2, set1));

}
// Jaccard similarity when one set is completely contained in the other
TEST(ScoresFixture, JaccardContained) {
    base::dynamic_bitset<> set1(5), set2(5);

    set1[0].set();
    set1[1].set();
    set1[2].set();
    set1[3].set();
    set2[0].set();
    set2[1].set();

    EXPECT_EQ(4, set1.count()) << "Set1 should have 4 members.";
    EXPECT_EQ(2, set2.count()) << "Set2 should have 2 members.";
    EXPECT_EQ(0.5, getJaccardSimilarity(set1, set2));
    EXPECT_EQ(0.5, getJaccardSimilarity(set2, set1));

}

// Jaccard similarity when one set is empty and the other not
TEST(ScoresFixture, JaccardOneEmpty) {
    base::dynamic_bitset<> set1(5), set2(5);

    set1[0].set();
    set1[1].set();
    set1[2].set();

    EXPECT_EQ(3, set1.count()) << "Set1 should have 3 members.";
    EXPECT_EQ(0, set2.count()) << "Set2 should have 0 members.";
    EXPECT_EQ(0.0, getJaccardSimilarity(set1, set2));
    EXPECT_EQ(0.0, getJaccardSimilarity(set2, set1));

}

// Jaccard similarity when they share a fraction of nodes
TEST(ScoresFixture, JaccardFraction) {
    base::dynamic_bitset<> set1(5), set2(5);

    set1[0].set();
    set1[1].set();
    set2[1].set();
    set2[2].set();

    EXPECT_EQ(2, set1.count()) << "Set1 should have 2 members.";
    EXPECT_EQ(2, set2.count()) << "Set2 should have 2 members.";
    EXPECT_NEAR(0.33, getJaccardSimilarity(set1, set2), 1e-2);
    EXPECT_NEAR(0.33, getJaccardSimilarity(set2, set1), 1e-2);

}
