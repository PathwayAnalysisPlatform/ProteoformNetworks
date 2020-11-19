#include "gtest/gtest.h"
#include "../Module.hpp"

TEST(ModuleSuite, DefaultConstructorThrowsErrorTest){
    ASSERT_THROW({Module module();}, std::runtime_error);
}


