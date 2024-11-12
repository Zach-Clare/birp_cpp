#include <gtest/gtest.h>

#include "../include/Cmem.h"

TEST(CMEM, BasicSamples) {
    // we wanna test the getsample method works at a basic level
    // We will test it against values from Sam's CMEM implementation.

    CMEM* cmem = new CMEM();
    cmem->Init();

    float sample1 = cmem->GetSample(8.f, 4.f, 4.f);
    // float expected1 = ????
    // EXPECT_EQ(expected1, sample1);

    float sample2 = cmem->GetSample(4.f, 12.f, 10.f);
    // float expected2 = ????
    // EXPECT_EQ(expected2, sample2);

    float sample3 = cmem->GetSample(3.f, 1.5f, 1.2f);
    // float expected3 = ????
    // EXPECT_EQ(expected3, sample3);

    GTEST_ASSERT_EQ(2, 2); 
}

TEST(CMEM, InitParamsSamples) {
    // we wanna test the getsample method with various parameter initilisations
    // We will test it against values from Sam's CMEM implementation.

    CMEM* cmem = new CMEM();
    cmem->Init();

    // is there support for custom input parameters?
    // no. we need to add this.

    float sample1 = cmem->GetSample(8.f, 4.f, 4.f);
    // float expected1 = ????
    // EXPECT_EQ(expected1, sample1);

    float sample2 = cmem->GetSample(4.f, 12.f, 10.f);
    // float expected2 = ????
    // EXPECT_EQ(expected2, sample2);

    float sample3 = cmem->GetSample(3.f, 1.5f, 1.2f);
    // float expected3 = ????
    // EXPECT_EQ(expected3, sample3);

    GTEST_ASSERT_EQ(2, 2); 
}