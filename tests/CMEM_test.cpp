#include <gtest/gtest.h>

#include "../include/Cmem.h"

TEST(CMEM, BasicSamples) {
    // we wanna test the getsample method works at a basic level
    // We will test it against values from Sam's CMEM implementation.

    CMEM* cmem = new CMEM();
    cmem->Init();

    float sample1 = cmem->GetSample(8.f, 4.f, 4.f);
    float expected1 = 0;
    EXPECT_EQ(expected1, sample1);

    float sample2 = cmem->GetSample(4.f, 12.f, 10.f);
    float expected2 = 3.05701769e-06;
    EXPECT_EQ(expected2, sample2);

    float sample3 = cmem->GetSample(3.f, 1.5f, 1.2f);
    float expected3 = 0;
    EXPECT_EQ(expected3, sample3);

    // checking for symmetry here
    float sample4 = cmem->GetSample(8.f, 10.f, 10.f);
    float expected4 = 8.15947317e-07;
    EXPECT_EQ(expected4, sample4);

    float sample5 = cmem->GetSample(8.f, -10.f, -10.f);
    float expected5 = 8.15947317e-07;
    EXPECT_EQ(expected5, sample5);

    float sample6 = cmem->GetSample(8.f, 10.f, -10.f);
    float expected6 = 8.15947317e-07;
    EXPECT_EQ(expected6, sample6);

    GTEST_ASSERT_EQ(2, 2); 
}

TEST(CMEM, InitParamsSamples) {
    // we wanna test the getsample method with various parameter initilisations

    std::vector<int> v = {400, 0, 0}; // 400, 0, 0
    std::vector<int> b = {0, 0, 4}; // 0, 0, 5
    float dipole = 0.f; // 0.f
    int p1 = 1; // 1
    int p2 = 4; // 3
    int p3 = 6; // 4
    float B = 3.f; // 2.f
    float alpha = 3.6f; // 2.5f
    float beta = -1.0f; // -1.6f
    

    CMEM* cmem = new CMEM();
    cmem->Init(true, v, b, dipole, p1, p2, p3, B, alpha, beta); // initialise with different paramets here and test against acceptable values below

    float sample1 = cmem->GetSample(8.f, 4.f, 4.f);
    float expected1 = 0;
    EXPECT_EQ(expected1, sample1);

    float sample2 = cmem->GetSample(4.f, 12.f, 10.f);
    float expected2 = 1.14234717e-06;
    EXPECT_EQ(expected2, sample2);

    float sample3 = cmem->GetSample(3.f, 1.5f, 1.2f);
    float expected3 = 0;
    EXPECT_EQ(expected3, sample3);

    // again check for symmetry
    // WARNING::::: These samples are identical to when we init with default parameters
    float sample4 = cmem->GetSample(8.f, 10.f, 10.f);
    float expected4 = 8.15947317e-07;
    EXPECT_EQ(expected4, sample4);

    float sample5 = cmem->GetSample(8.f, -10.f, -10.f);
    float expected5 = 8.15947317e-07;
    EXPECT_EQ(expected5, sample5);

    float sample6 = cmem->GetSample(8.f, 10.f, -10.f);
    float expected6 = 8.15947317e-07;
    EXPECT_EQ(expected6, sample6);
}