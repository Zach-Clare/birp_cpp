/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


#include <gtest/gtest.h>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>

#include "../include/Camera.h"
#include <Cmem.h>

int main()
{
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}

class CubeFake : public Space {
    public:
        inline void Init() { };
        inline float GetSample(float x, float y, float z) { return 1.f; }
};

TEST(Camera, ToDatCreatesFiles) {
    Space* cube = new CubeFake;
    float pixel_size_deg = 0.25f;
    float plot_fov_h = 36.f;
    float plot_fov_w = 36.f;
    Camera camera(*cube, pixel_size_deg, plot_fov_h, plot_fov_w);
    float x = 1.f;
    float y = 1.f;
    float z = 1.f;
    camera.SetPosition(x, y, z);
    camera.SetAim(1, 0, 0);
    camera.Render();
    std::filesystem::remove("/home/zc/code/birp/cpp/data/output/testout.dat");

    camera.ToDat("/home/zc/code/birp/cpp/data/output/testout");

    FILE* f = fopen("/home/zc/code/birp/cpp/data/output/testout.dat", "r");
    ASSERT_TRUE(f != NULL);
    fclose(f);
}

TEST(Camera, ToFITSCreatesFile) {
    // We currently aren't testing if the FITS file is malformed or otherwise wrong.
    // We are just testing if it creates a file with a fits extension in the specified place.
    Space* cube = new CubeFake;
    float pixel_size_deg = 0.25f;
    float plot_fov_h = 36.f;
    float plot_fov_w = 36.f;
    Camera camera(*cube, pixel_size_deg, plot_fov_h, plot_fov_w);
    float x = 1.f;
    float y = 1.f;
    float z = 1.f;
    camera.SetPosition(x, y, z);
    camera.SetAim(1, 0, 0);
    camera.Render();
    std::filesystem::remove("/home/zc/code/birp/cpp/data/output/testout.fits");

    camera.ToFITS("/home/zc/code/birp/cpp/data/output/testout");

    FILE* f = fopen("/home/zc/code/birp/cpp/data/output/testout.fits", "r");
    ASSERT_TRUE(f != NULL);
    fclose(f);
}

////////
// I don't like the idea of taking a full render and checking individual pixels as the basis of a test
// We're trying to check that an image comes out? Or that CMEM is correctly implemented? Or that a rendering error hasn't botched the values?
////////
// TEST(Camera, RendersAccurateMHD) {
//     DataCube* cube = new DataCube;
//     Space* space;
//     space = cube;
//     std::string path = "/home/zc/code/birp/cpp/data/px_uni_1000.dat";
//     cube->Init(path, false);
//     cube->SetTrilinear(false);
//     float pixel_size_deg = 0.25f;
//     float plot_fov_h = 36.f;
//     float plot_fov_w = 36.f;
//     Camera camera(*cube, pixel_size_deg, plot_fov_h, plot_fov_w);
//     float x = 1.f;
//     float y = 1.f;
//     float z = 1.f;
//     camera.SetPosition(x, y, z);
//     camera.SetAim(1, 0, 0);
//     camera.Render();
//     std::filesystem::remove("/home/zc/code/birp/cpp/data/output/testoutMHD.dat"); // remove residual files

//     camera.ToDat("/home/zc/code/birp/cpp/data/output/testoutMHD"); // out in readable format
        
//     std::fstream fin; 
//     std::fstream infile("/home/zc/code/birp/cpp/data/output/testoutMHD.dat"); // read output from birp

//     std::vector<std::string> row;
//     std::string temp;

//     while (std::getline(infile, temp, ',')) {
//         row.push_back(temp); // loop through and save all values
//     }

//     // pluck a series of random values to test the render was accurate
//     float sample1 = std::stof(row[93]);
//     float expected1 = 30.806702f;
//     EXPECT_EQ(expected1, sample1);

//     float sample2 = std::stof(row[751]);
//     float expected2 = 24.165386f;
//     EXPECT_EQ(expected2, sample2);

//     // EXPECT_EQ(2, 2);
// }

// TEST(Camera, RendersAccurateCMEM) {
//     CMEM* cube = new CMEM;
//     Space* space;
//     space = cube;
//     cube->Init();
//     float pixel_size_deg = 0.25f;
//     float plot_fov_h = 36.f;
//     float plot_fov_w = 36.f;
//     Camera camera(*cube, pixel_size_deg, plot_fov_h, plot_fov_w);
//     float x = 1.f;
//     float y = 1.f;
//     float z = 1.f;
//     camera.SetPosition(x, y, z);
//     camera.SetAim(1, 0, 0);
//     camera.Render();
//     std::filesystem::remove("/home/zc/code/birp/cpp/data/output/testoutCMEM.dat"); // remove residual files

//     camera.ToDat("/home/zc/code/birp/cpp/data/output/testoutCMEM"); // out in readable format
        
//     std::fstream fin; 
//     std::fstream infile("/home/zc/code/birp/cpp/data/output/testoutCMEM.dat"); // read output from birp

//     std::vector<std::string> row;
//     std::string temp;

//     while (std::getline(infile, temp, ',')) {
//         row.push_back(temp); // loop through and save all values
//     }

//     // pluck a series of random values to test the render was accurate
//     float sample1 = std::stof(row[93]);
//     float expected1 = 0.696530998f;
//     EXPECT_EQ(expected1, sample1);

//     float sample2 = std::stof(row[751]);
//     float expected2 = 1.03772402f;
//     EXPECT_EQ(expected2, sample2);
// }