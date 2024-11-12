#include <gtest/gtest.h>

#include "../include/Camera.h"

TEST(Camera, init) {
    GTEST_ASSERT_EQ(2, 2);
}

int main()
{
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}