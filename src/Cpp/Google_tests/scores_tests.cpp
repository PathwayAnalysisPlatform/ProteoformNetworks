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

// Overlap similarity: both sets empty
TEST(ScoresFixtyre, OverlapBothEmpty) {
    base::dynamic_bitset<> set1(5), set2(5);

    EXPECT_EQ(0, set1.count()) << "Set1 should be empty.";
    EXPECT_EQ(0, set2.count()) << "Set2 should be empty.";
    EXPECT_EQ(1, getOverlapSimilarity(set1, set2)) << "Both sets are empty then overlap should be 1.";
    EXPECT_EQ(1, getOverlapSimilarity(set2, set1)) << "Both sets are empty then overlap should be 1.";
}

// Overlap similarity: disjoint sets
TEST(ScoresFixture, OverlapDisjoint) {
    base::dynamic_bitset<> set1(5), set2(5);

    set1[0].set();
    set1[1].set();
    set2[2].set();
    set2[3].set();

    EXPECT_EQ(2, set1.count()) << "Set1 should have 2 elements.";
    EXPECT_EQ(2, set2.count()) << "Set2 should have 2 elements.";
    EXPECT_EQ(0, getOverlapSimilarity(set1, set2)) << "Sets should not overlap: then score 0.";
    EXPECT_EQ(0, getOverlapSimilarity(set2, set1)) << "Sets should not overlap: then score 0.";
}

// Overlap similarity: both same set
TEST(ScoresFixture, BothSameSet){
    base::dynamic_bitset<> set1(5), set2(5);
    set1[1].set();
    set1[3].set();
    set2[1].set();
    set2[3].set();

    EXPECT_EQ(2, set1.count()) << "Set1 should have 2 elements.";
    EXPECT_EQ(2, set2.count()) << "Set2 should have 2 elements.";
    EXPECT_EQ(1.0, getOverlapSimilarity(set1, set2)) << "Sets are the same, then should score 1.";
    EXPECT_EQ(1.0, getOverlapSimilarity(set2, set1)) << "Sets are the same, then should score 1.";
}

// Overlap similarity: one contained in the other set
TEST(ScoresFixture, OnIsSubset) {
    base::dynamic_bitset<> set1(5), set2(5);
    set1[0].set();
    set1[1].set();
    for(int I = 0; I < 5; I++){
        set2[I].set();
    }

    EXPECT_EQ(2, set1.count()) << "Set1 should have 2 elements.";
    EXPECT_EQ(5, set2.count()) << "Set2 should have 5 elements.";
    EXPECT_EQ(1.0, getOverlapSimilarity(set1, set2)) << "The smaller set is contained in the other set, then score 1.";
    EXPECT_EQ(1.0, getOverlapSimilarity(set2, set1)) << "The smaller set is contained in the other set, then score 1.";
}
// Overlap similarity: one set is empty the other is not
TEST(ScoresFixture, OneEmpty) {
    base::dynamic_bitset<> set1(5), set2(5);
    set1[1].set();
    set1[2].set();
    set1[3].set();

    EXPECT_EQ(3, set1.count()) << "Set1 should have 3 elements.";
    EXPECT_EQ(0, set2.count()) << "Set2 should be empty.";
    EXPECT_EQ(1, getOverlapSimilarity(set1, set2)) << "One set is empty, then score 0.";
    EXPECT_EQ(1, getOverlapSimilarity(set2, set1)) << "One set is empty, then score 0.";
}

// Overlap similarity: sets share fraction of elements
TEST(ScoresFixture, ShareFraction){
    base::dynamic_bitset<> set1(5), set2(5);
    for(int I = 0; I < 4; I++)  // 0, 1, 2, 3
        set1[I].set();
    for(int I = 1; I < 5; I++)  // 1, 2, 3, 4
        set2[I].set();

    EXPECT_EQ(4, set1.count()) << "Set1 should have 4 elements.";
    EXPECT_EQ(4, set2.count()) << "Set2 should have 4 elements.";
    EXPECT_NEAR(0.75, getOverlapSimilarity(set1, set2), 1e-2) << "Sets should share three quarters of members.";
    EXPECT_NEAR(0.75, getOverlapSimilarity(set2, set1), 1e-2) << "Sets should share three quarters of members.";
}

// Overlap similarity: sets share all the elements of the smaller set
TEST(ScoresFixture, ShareSmallSet){
    base::dynamic_bitset<> set1(5), set2(5);
    for(int I = 0; I <= 2; I++)  // 0, 1, 2
        set1[I].set();
    for(int I = 0; I < 5; I++)  // 1, 2, 3, 4, 5
        set2[I].set();

    EXPECT_EQ(3, set1.count()) << "Set1 should have 3 elements.";
    EXPECT_EQ(5, set2.count()) << "Set2 should have 5 elements.";
    EXPECT_NEAR(1, getOverlapSimilarity(set1, set2), 1e-2) << "Sets should share all smaller set members.";
    EXPECT_NEAR(1, getOverlapSimilarity(set2, set1), 1e-2) << "Sets should share all smaller set members.";
}
