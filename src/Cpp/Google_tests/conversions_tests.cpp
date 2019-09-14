#include <conversions.hpp>
#include "gtest/gtest.h"
#include "types.hpp"

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
    for(const auto &vec_element : myStrVec){
        EXPECT_TRUE(myStrSet.find(vec_element) != myStrSet.end());
    }
}
