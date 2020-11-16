#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../Interactome.hpp"

class ModuleCreatorFixture : public ::testing::Test {

protected:
    virtual void SetUp() override {

    }

    virtual void TearDown() override {

    }

    std::string_view file_phegeni;
    Interactome interactome;
    const std::string output_path;
};

TEST(CreateModulesSuite, CreateGeneModuleTest){
    ASSERT_TRUE(false);
}

TEST(CreateModulesSuite, CreateProteinModuleTest){
    ASSERT_TRUE(false);
}

TEST(CreateModulesSuite, CreateProetoformModuleTest){
    ASSERT_TRUE(false);
}

TEST(CreateModulesSuite, SaveModuleTest){
    ASSERT_TRUE(false);
}