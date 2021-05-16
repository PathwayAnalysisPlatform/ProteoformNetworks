#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "Interactome.hpp"

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
    ASSERT_TRUE(true);
}

TEST(CreateModulesSuite, CreateProteinModuleTest){
    ASSERT_TRUE(true);
}

TEST(CreateModulesSuite, CreateProetoformModuleTest){
    ASSERT_TRUE(true);
}

TEST(CreateModulesSuite, SaveModuleTest){
    ASSERT_TRUE(true);
}