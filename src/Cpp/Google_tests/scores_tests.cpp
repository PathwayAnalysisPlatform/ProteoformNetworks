#include <bitset>
#include <algorithm>
#include "gtest/gtest.h"
#include "scores.hpp"

class ScoresFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        // Create 4 sets of 5 members.
        sets["A"] = base::dynamic_bitset<>(5);  // 0, 1, 2, 3, 4
        sets["B"] = base::dynamic_bitset<>(5);  // 0, 1, 2
        sets["C"] = base::dynamic_bitset<>(5);  // 3, 4
        sets["D"] = base::dynamic_bitset<>(5);  // 1, 2, 3

        for (int I = 0; I <= 4; I++)
            sets["A"][I].set();
        for (int I = 0; I <= 2; I++)
            sets["B"][I].set();
        for (int I = 3; I <= 4; I++)
            sets["C"][I].set();
        for (int I = 1; I <= 3; I++)
            sets["D"][I].set();
    }

    msb sets;
};

// Jaccard similarity when both sets are empty
TEST_F(ScoresFixture, JaccardBothEmpty) {
    base::dynamic_bitset<> set1, set2;

    EXPECT_EQ(0, set1.count()) << "The set1 should be empty.";
    EXPECT_EQ(0, set2.count()) << "The set2 should be empty.";
    EXPECT_EQ(1, getJaccardSimilarity(set1, set2));
    EXPECT_EQ(1, getJaccardSimilarity(set2, set1));

}

// Jaccard similarity when the sets are disjoint (completely different)
TEST_F(ScoresFixture, JaccardDisjoint) {
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
TEST_F(ScoresFixture, JaccardBothSame) {
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
TEST_F(ScoresFixture, JaccardContained) {
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
TEST_F(ScoresFixture, JaccardOneEmpty) {
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
TEST_F(ScoresFixture, JaccardFraction) {
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
TEST_F(ScoresFixture, OverlapBothEmpty) {
    base::dynamic_bitset<> set1(5), set2(5);

    EXPECT_EQ(0, set1.count()) << "Set1 should be empty.";
    EXPECT_EQ(0, set2.count()) << "Set2 should be empty.";
    EXPECT_EQ(1, getOverlapSimilarity(set1, set2)) << "Both sets are empty then overlap should be 1.";
    EXPECT_EQ(1, getOverlapSimilarity(set2, set1)) << "Both sets are empty then overlap should be 1.";
}

// Overlap similarity: disjoint sets
TEST_F(ScoresFixture, OverlapDisjoint) {
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
TEST_F(ScoresFixture, BothSameSet) {
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
TEST_F(ScoresFixture, OnIsSubset) {
    base::dynamic_bitset<> set1(5), set2(5);
    set1[0].set();
    set1[1].set();
    for (int I = 0; I < 5; I++) {
        set2[I].set();
    }

    EXPECT_EQ(2, set1.count()) << "Set1 should have 2 elements.";
    EXPECT_EQ(5, set2.count()) << "Set2 should have 5 elements.";
    EXPECT_EQ(1.0, getOverlapSimilarity(set1, set2)) << "The smaller set is contained in the other set, then score 1.";
    EXPECT_EQ(1.0, getOverlapSimilarity(set2, set1)) << "The smaller set is contained in the other set, then score 1.";
}
// Overlap similarity: one set is empty the other is not
TEST_F(ScoresFixture, OneEmpty) {
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
TEST_F(ScoresFixture, ShareFraction) {
    base::dynamic_bitset<> set1(5), set2(5);
    for (int I = 0; I < 4; I++)  // 0, 1, 2, 3
        set1[I].set();
    for (int I = 1; I < 5; I++)  // 1, 2, 3, 4
        set2[I].set();

    EXPECT_EQ(4, set1.count()) << "Set1 should have 4 elements.";
    EXPECT_EQ(4, set2.count()) << "Set2 should have 4 elements.";
    EXPECT_NEAR(0.75, getOverlapSimilarity(set1, set2), 1e-2) << "Sets should share three quarters of members.";
    EXPECT_NEAR(0.75, getOverlapSimilarity(set2, set1), 1e-2) << "Sets should share three quarters of members.";
}

// Overlap similarity: sets share all the elements of the smaller set
TEST_F(ScoresFixture, ShareSmallSet) {
    base::dynamic_bitset<> set1(5), set2(5);
    for (int I = 0; I <= 2; I++)  // 0, 1, 2
        set1[I].set();
    for (int I = 0; I < 5; I++)  // 1, 2, 3, 4, 5
        set2[I].set();

    EXPECT_EQ(3, set1.count()) << "Set1 should have 3 elements.";
    EXPECT_EQ(5, set2.count()) << "Set2 should have 5 elements.";
    EXPECT_NEAR(1, getOverlapSimilarity(set1, set2), 1e-2) << "Sets should share all smaller set members.";
    EXPECT_NEAR(1, getOverlapSimilarity(set2, set1), 1e-2) << "Sets should share all smaller set members.";
}

// Get scores for pairs of sets: Correct number of pairs
TEST_F(ScoresFixture, AllPairsCorrectPairs) {
    // Create 4 sets of 5 members.
    msb sets;
    sets["A"] = base::dynamic_bitset<>(5);
    sets["B"] = base::dynamic_bitset<>(5);
    sets["C"] = base::dynamic_bitset<>(5);
    sets["D"] = base::dynamic_bitset<>(5);

    auto scores = getScores(sets, getOverlapSimilarity);

    EXPECT_EQ((sets.size() * (sets.size() - 1)) / 2, scores.size());
}

// Get scores for pairs of sets: Correct values
TEST_F(ScoresFixture, AllPairsCorrectValues) {
    um<std::string, double> scores = getScores(sets, getOverlapSimilarity);

    for (const auto &score_entry : scores) {
        std::cout << "Score " << score_entry.first << ": " << score_entry.second << std::endl;
    }

    // 0: (A, B)
    // 1: (A, C)
    // 2: (A, D)
    // 3: (B, C)
    // 4: (B, D)
    // 5: (C, D)

    EXPECT_EQ(1, scores["A\tB"]);    // A B
    EXPECT_EQ(0, scores["B\tC"]);    // B C
    EXPECT_NEAR(0.66, scores["B\tD"], 1e-1); // B D
    EXPECT_EQ(0.5, scores["C\tD"]);  // C D
}


TEST_F(ScoresFixture, WriteScores) {
    // Create example file
    std::string file_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    file_name += "_scores.txt";
    writeScores(sets, getScores(sets, getOverlapSimilarity), file_name);

    // Open file
    std::ifstream scores_file(file_name);
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(scores_file, line)) {
        lines.push_back(line);
    }


    EXPECT_EQ(7, lines.size()); // Check number of lines is correct
    EXPECT_EQ("TRAIT1\tTRAIT2\tSCORE", lines[0]); // Check header is there
    // Check every row has 3 columns
    for (int I = 0; I < lines.size(); I++) {
        size_t appeareances = std::count(lines[I].begin(), lines[I].end(), '\t');
        EXPECT_EQ(2, (int) appeareances);
    }
    // Check values are correct
    EXPECT_EQ("A\tB\t1", lines[1]);
    EXPECT_EQ("B\tC\t0", lines[4]);
    EXPECT_EQ("B\tD\t0.66667", lines[5]);
    EXPECT_EQ("C\tD\t0.5", lines[6]);

    // Delete example file
    int n = file_name.length();
    char file_name_char[n + 1];
    strcpy(file_name_char, file_name.c_str());
    std::remove(file_name_char);
}
