/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


#include <gtest/gtest.h>
#include <string>

#include "DataCube.h"

TEST(DataCube, correctSample) {
    DataCube* cube = new DataCube;
    std::string path = "../data/px_uni_1000.dat";
    cube->Init(path, 82, false);
    cube->SetTrilinear(false);

    float sample = cube->GetSample(9.f, 0.f, 0.f);
    float expected = 0.000144541002;
    EXPECT_EQ(expected, sample);

    float sample2 = cube->GetSample(4.f, 9.4f, 6.f);
    float expected2 = 3.21660991e-05;
    EXPECT_EQ(expected2, sample2);
}

TEST(DataCube, correctSampleTrilinear) {
    DataCube* cube = new DataCube;
    std::string path = "../data/px_uni_1000.dat";
    cube->Init(path, 82, false);
    cube->SetTrilinear(true);

    float sample = cube->GetSample(9.f, 0.f, 0.f);
    float expected = 0.000144541002;
    EXPECT_EQ(expected, sample);

    float sample2 = cube->GetSample(4.f, 9.4f, 6.f);
    float expected2 = 3.11356634e-05;
    EXPECT_EQ(expected2, sample2);
}

TEST(DataCube, readDataFile) {
    DataCube* cube = new DataCube;
    std::string path = "../data/px_uni_1000.dat";
    cube->Init(path, 82, false);

    EXPECT_EQ(101, cube->size.x);
    EXPECT_EQ(161, cube->size.y);
    EXPECT_EQ(161, cube->size.z);
}
