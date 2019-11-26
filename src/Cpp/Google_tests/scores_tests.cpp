#include <bitset>
#include <algorithm>
#include <utility>
#include "gtest/gtest.h"
#include "scores.hpp"

class ScoresFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        groups = createBimap({"A", "B", "C", "D"});

        // Create 4 sets of 5 members.
        sets.push_back(base::dynamic_bitset<>(5));  // 0, 1, 2, 3, 4
        sets.push_back(base::dynamic_bitset<>(5));  // 0, 1, 2
        sets.push_back(base::dynamic_bitset<>(5));  // 3, 4
        sets.push_back(base::dynamic_bitset<>(5));  // 1, 2, 3

        for (int I = 0; I <= 4; I++)
            sets[groups.str_to_int["A"]][I].set();
        for (int I = 0; I <= 2; I++)
            sets[groups.str_to_int["B"]][I].set();
        for (int I = 3; I <= 4; I++)
            sets[groups.str_to_int["C"]][I].set();
        for (int I = 1; I <= 3; I++)
            sets[groups.str_to_int["D"]][I].set();

        the_modules.group_to_members = sets;

    }

    vb sets;
    modules the_modules;
    bimap_str_int groups, members;
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
    vb sets;
    sets.push_back(base::dynamic_bitset<>(5));
    sets.push_back(base::dynamic_bitset<>(5));
    sets.push_back(base::dynamic_bitset<>(5));
    sets.push_back(base::dynamic_bitset<>(5));

    auto scores = getScores(sets, getOverlapSimilarity, 0, 6);

    EXPECT_EQ((sets.size() * (sets.size() - 1)) / 2, scores.size());
}

// Get scores for pairs of sets: Correct values
TEST_F(ScoresFixture, AllPairsCorrectValues) {
    pair_map<double> scores = getScores(sets, getOverlapSimilarity, 0, 6);

    for (const auto &score_entry : scores) {
        std::cout << "Score " << groups.int_to_str[score_entry.first.first] << '\t'
                  << groups.int_to_str[score_entry.first.second] << ": " << score_entry.second << std::endl;
    }

    // 0: (A, B)
    // 1: (A, C)
    // 2: (A, D)
    // 3: (B, C)
    // 4: (B, D)
    // 5: (C, D)

    EXPECT_EQ(1, scores[std::make_pair(groups.str_to_int["A"], groups.str_to_int["B"])]);    // A B
    EXPECT_EQ(0, scores[std::make_pair(groups.str_to_int["B"], groups.str_to_int["C"])]);    // B C
    EXPECT_NEAR(0.66, scores[std::make_pair(groups.str_to_int["B"], groups.str_to_int["D"])], 1e-1); // B D
    EXPECT_EQ(0.5, scores[std::make_pair(groups.str_to_int["C"], groups.str_to_int["D"])]);  // C D
}


TEST_F(ScoresFixture, WriteScores) {
    // Create example file
    std::string file_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    file_name += "_scores.txt";

    std::vector<std::string> features_labels = {"SCORE"};
    std::vector<pair_map<double>> features = {getScores(sets, getOverlapSimilarity, 0, 10)};
    writeScores(groups, the_modules, features_labels, features, file_name);

    // Open file
    std::ifstream scores_file(file_name);
    std::string line;
    std::set<std::string> lines;
    while (std::getline(scores_file, line)) {
        lines.insert(line);
    }

    EXPECT_EQ(6, lines.size()); // Check number of lines is correct
    EXPECT_TRUE(lines.find("TRAIT1\tTRAIT2\tSCORE") != lines.end()); // Check header is there
    // Check every row has 3 columns
    for(const auto& line : lines){
        size_t appeareances = std::count(line.begin(), line.end(), '\t');
        EXPECT_EQ(2, (int) appeareances);
    }
    // Check values are correct
    EXPECT_TRUE(lines.find("A\tB\t1") != lines.end());
    EXPECT_TRUE(lines.find("B\tC\t0") == lines.end());
    EXPECT_TRUE(lines.find("B\tD\t0.66667") != lines.end());
    EXPECT_TRUE(lines.find("C\tD\t0.5") != lines.end());

    // Delete example file
    int n = file_name.length();
    char file_name_char[n + 1];
    strcpy(file_name_char, file_name.c_str());
    std::remove(file_name_char);
}

/*
 * Fixture to test three methods:
 * - Interface size calculation in nodes
 * - Interface size calculation in edges
 * - Interface size calculation for all pairs of sets
 */
class InterfaceSizeFixture : public ::testing::Test {

protected:
    virtual void SetUp() {
        groups = createBimap({"A", "B", "C", "D"});
//
//        // Create 4 sets of 5 members.
//        sets.push_back(base::dynamic_bitset<>(5));  // 0, 1, 2, 3, 4
//        sets.push_back(base::dynamic_bitset<>(5));  // 0, 1, 2
//        sets.push_back(base::dynamic_bitset<>(5));  // 3, 4
//        sets.push_back(base::dynamic_bitset<>(5));  // 1, 2, 3
//
//        for (int I = 0; I <= 4; I++)
//            sets[groups.str_to_int["A"]][I].set();
//        for (int I = 0; I <= 2; I++)
//            sets[groups.str_to_int["B"]][I].set();
//        for (int I = 3; I <= 4; I++)
//            sets[groups.str_to_int["C"]][I].set();
//        for (int I = 1; I <= 3; I++)
//            sets[groups.str_to_int["D"]][I].set();
//
//        the_modules.group_to_members = sets;

    }

//
//    vb sets;
//    modules the_modules;
    bimap_str_int groups, members;
};

TEST_F(InterfaceSizeFixture, SetsDoNotOverlapAndNotConnectedSoNoInterface) {
    // Create sets: vertices and edges
    base::dynamic_bitset<> V1 = base::dynamic_bitset<>(5);
    base::dynamic_bitset<> V2 = base::dynamic_bitset<>(5);
    vusi E(5, std::unordered_set<int>());

    V1[0].set();
    V1[1].set();
    V2[2].set();
    V2[3].set();
    V2[4].set();

    E[0].insert(1);
    E[0].insert(0);
    E[2].insert({3, 4});
    E[3].insert({2, 4});
    E[4].insert({2, 3});

    int size_nodes = calculate_interface_size_nodes(V1, V2, E);
    int size_edges = calculate_interface_size_edges(V1, V2, E);
    EXPECT_EQ(0, size_nodes);
    EXPECT_EQ(0, size_edges);
}

TEST_F(InterfaceSizeFixture, OneSetIsEmptySoNoInterface) {
    base::dynamic_bitset<> V1 = base::dynamic_bitset<>(5);
    base::dynamic_bitset<> V2 = base::dynamic_bitset<>(5);
    vusi E(5, std::unordered_set<int>());

    V2[2].set();
    V2[3].set();
    V2[4].set();

    E[2].insert({3, 4});
    E[3].insert({2, 4});
    E[4].insert({2, 3});

    int size_nodes = calculate_interface_size_nodes(V1, V2, E);
    int size_edges = calculate_interface_size_edges(V1, V2, E);
    EXPECT_EQ(0, size_nodes);
    EXPECT_EQ(0, size_edges);
}

TEST_F(InterfaceSizeFixture, BothSetsAreEmptySoNoInterface) {
    base::dynamic_bitset<> V1 = base::dynamic_bitset<>(5);
    base::dynamic_bitset<> V2 = base::dynamic_bitset<>(5);
    vusi E(5, std::unordered_set<int>());

    int size_nodes = calculate_interface_size_nodes(V1, V2, E);
    int size_edges = calculate_interface_size_edges(V1, V2, E);
    EXPECT_EQ(0, size_nodes);
    EXPECT_EQ(0, size_edges);
}

TEST_F(InterfaceSizeFixture, SetsAreConnectedAndNotOverlap) {
    base::dynamic_bitset<> V1 = base::dynamic_bitset<>(8);
    base::dynamic_bitset<> V2 = base::dynamic_bitset<>(8);
    vusi E(8, std::unordered_set<int>());

    for (int I : {0, 1, 2, 3})
        V1[I].set();
    for (int I : {4, 5, 6, 7})
        V2[I].set();

    E[0].insert({1, 2});
    E[1].insert({0, 3, 4});
    E[2].insert({2, 3});
    E[3].insert({1, 2, 4, 6});
    E[4].insert({1, 3, 5, 6});
    E[5].insert({4, 7});
    E[6].insert({3, 4, 7});
    E[7].insert({5, 6});

    int size_nodes = calculate_interface_size_nodes(V1, V2, E);
    int size_edges = calculate_interface_size_edges(V1, V2, E);
    EXPECT_EQ(4, size_nodes);
    EXPECT_EQ(3, size_edges);
}

TEST_F(InterfaceSizeFixture, SetsOverlapAtOneVertexAndConnectedOnlyFromOverlapVertices) {
    base::dynamic_bitset<> V1 = base::dynamic_bitset<>(7);
    base::dynamic_bitset<> V2 = base::dynamic_bitset<>(7);
    vusi E(7, std::unordered_set<int>());

    for (int I : {0, 1, 2, 3})
        V1[I].set();
    for (int I : {3, 4, 5, 6})
        V2[I].set();

    E[0].insert({1, 2});
    E[1].insert({0, 3});
    E[2].insert({0, 3});
    E[3].insert({1, 2, 4, 5});
    E[4].insert({3, 5, 6});
    E[5].insert({4, 6});
    E[6].insert({4, 5});

    int size_nodes = calculate_interface_size_nodes(V1, V2, E);
    int size_edges = calculate_interface_size_edges(V1, V2, E);
    EXPECT_EQ(5, size_nodes);
    EXPECT_EQ(4, size_edges);
}

TEST_F(InterfaceSizeFixture, SetsOverlapAtMultipleVerticesAndConnectOnlyFromOverlapVertices) {
    base::dynamic_bitset<> V1 = base::dynamic_bitset<>(14);
    base::dynamic_bitset<> V2 = base::dynamic_bitset<>(14);
    vusi E(14, std::unordered_set<int>());

    for (int I : {1, 2, 3, 4, 5, 6, 7})
        V1[I].set();
    for (int I : {0, 2, 4, 7, 8, 9, 10, 11, 12, 13})
        V2[I].set();

    E[0].insert({8, 9, 10});
    E[1].insert({2, 5});
    E[2].insert({1, 3, 7});
    E[3].insert({2, 5, 6});
    E[4].insert({6, 7, 8});
    E[5].insert({1, 3});
    E[6].insert({3, 4});
    E[7].insert({2, 4, 10});
    E[8].insert({0, 4, 9});
    E[9].insert({0, 8});
    E[10].insert({0, 7});
    E[11].insert({12, 13});
    E[12].insert({11, 13});
    E[13].insert({11, 12});

    int size_nodes = calculate_interface_size_nodes(V1, V2, E);
    int size_edges = calculate_interface_size_edges(V1, V2, E);
    EXPECT_EQ(8, size_nodes);
    EXPECT_EQ(7, size_edges);
}

TEST_F(InterfaceSizeFixture, SetsOverlapAndConnectedFromOverlapAndNonOverlappingVertices) {
    base::dynamic_bitset<> V1 = base::dynamic_bitset<>(14);
    base::dynamic_bitset<> V2 = base::dynamic_bitset<>(14);
    vusi E(14, std::unordered_set<int>());

    for (int I : {1, 2, 3, 4, 5, 6, 7})
        V1[I].set();
    for (int I : {0, 2, 4, 7, 8, 9, 10, 11, 12, 13})
        V2[I].set();

    E[0].insert({1, 8, 9, 10});
    E[1].insert({0, 2, 5});
    E[2].insert({1, 3, 7});
    E[3].insert({2, 5, 6});
    E[4].insert({6, 7, 8});
    E[5].insert({1, 3});
    E[6].insert({3, 4, 9});
    E[7].insert({2, 4, 10});
    E[8].insert({0, 4, 9});
    E[9].insert({0, 6, 8});
    E[10].insert({0, 7});
    E[11].insert({12, 13});
    E[12].insert({11, 13});
    E[13].insert({11, 12});

    int size_nodes = calculate_interface_size_nodes(V1, V2, E);
    int size_edges = calculate_interface_size_edges(V1, V2, E);
    EXPECT_EQ(10, size_nodes);
    EXPECT_EQ(9, size_edges);
}
