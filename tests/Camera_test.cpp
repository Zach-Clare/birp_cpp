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
#include <algorithm>

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
    std::filesystem::remove("../data/output/testout.dat");

    camera.ToDat("../data/output/testout");

    FILE* f = fopen("../data/output/testout.dat", "r");
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
    std::filesystem::remove("../data/output/testout.fits");

    camera.ToFITS("../data/output/testout");

    FILE* f = fopen("../data/output/testout.fits", "r");
    ASSERT_TRUE(f != NULL);
    fclose(f);
}

////////
// I don't like the idea of taking a full render and checking individual pixels as the basis of a test
// We're trying to check that an image comes out? Or that CMEM is correctly implemented? Or that a rendering error hasn't botched the values?
////////
// Okay so I used to believe that comment above until an issue calculating the ingress and index in dataCube::GetSample() caused an issue with a researcher's work 
// and I'm no longer a purist, I have graduated into a realist. Oh no, entire renders! Improperly siloed unit tests! The unit tests are too useful! The horror!!!
// I'm being a bit facetious, I understand why specific unit tests are better and why this here isn't exactly a unit tests when everything else isn't mocked.
// Actually...why am I not mocking this? That's a great question. The reason I'm bringing back this test is to prevent a weird dragon-shaped magnetopause being rendered
// due to an issue with the cube index sampling method. I can test that sampling method with a bunch of edge cases and all that, but I think it's far better to just use
// some real data and make sure the same image is rendered every time. I'm not sure how much i agree with all this, but at some point we have to stop debating and do
// something that at least works. So here it is. 
////////
TEST(Camera, RendersAccurateMHD) {
    DataCube* cube = new DataCube;
    Space* space;
    space = cube;
    std::string path = "../data/px_uni_1000.dat";
    cube->Init(path, false);
    cube->SetTrilinear(false);
    float pixel_size_deg = 0.25f;
    float plot_fov_h = 36.f;
    float plot_fov_w = 36.f;
    Camera camera(*cube, pixel_size_deg, plot_fov_h, plot_fov_w);
    float x = 5.f;
    float y = 20.f;
    float z = 0.f;
    camera.SetPosition(x, y, z);
    camera.SetAim(6, 0, 0);
    camera.Render();
    std::filesystem::remove("../data/output/testoutMHD.dat"); // remove residual files

    camera.ToDat("../data/output/testoutMHD"); // out in readable format
    camera.ToFITS("../data/output/testoutMHD"); // out in readable format
        
    std::fstream fin; 
    std::fstream output_file("../data/output/testoutMHD.dat"); // read output from birp
    std::fstream truth_file("../data/test/test_MHD_static.dat"); // read output from a known good render
    std::fstream incorrect_file("../data/test/test_MHD_static_incorrect.dat"); // just to be sure, try a known bad render

    std::vector<std::string> output_image;
    std::vector<std::string> truth_image;
    std::vector<std::string> incorrect_image;
    std::string output_row;
    std::string truth_row;
    std::string incorrect_row;

    while (std::getline(output_file, output_row, ',')) {
        output_image.push_back(output_row);
    }
    while (std::getline(truth_file, truth_row, ',')) {
        truth_image.push_back(truth_row);
    }
    while (std::getline(incorrect_file, incorrect_row, ',')) {
        incorrect_image.push_back(incorrect_row);
    }

    EXPECT_TRUE(std::equal( truth_image.begin(), truth_image.end(), output_image.begin() ));
    EXPECT_FALSE(std::equal( incorrect_image.begin(), incorrect_image.end(), output_image.begin() ));
}

TEST(Camera, RendersAccurateCMEM) {
    CMEM* cube = new CMEM;
    Space* space;
    space = cube;
    //CMEM init values
    std::vector<int> v = {400, 0, 0};
    std::vector<int> b = {-5, 0, 0};
    float dipole = 30.f;
    float p0 = 0.67f;
    int p1 = 1;
    int p2 = 3;
    int p3 = 4;
    float B = 2.f;
    float alpha = 2.5f;
    float beta = -1.6f;
    float bs = 12.6400003f;
    float A1 = 7.2e-06f;
    float A2 = 3.5000000000000004e-06f;
    float ay_bs = NULL;
    float az_bs = NULL;
    float density = NULL;
    cube->Init(false, v, b, dipole, p0, p1, p2, p3, B, alpha, beta, bs, A1, A2, ay_bs, az_bs, density);
    float pixel_size_deg = 0.25f;
    float plot_fov_h = 36.f;
    float plot_fov_w = 36.f;
    Camera camera(*cube, pixel_size_deg, plot_fov_h, plot_fov_w);
    float x = 5.f;
    float y = 20.f;
    float z = 0.f;
    camera.SetPosition(x, y, z);
    camera.SetAim(6, 0, 0);
    camera.Render();
    std::filesystem::remove("../data/output/testoutCMEM.dat"); // remove residual files

    camera.ToDat("../data/output/testoutCMEM"); // out in readable format
    camera.ToFITS("../data/output/testoutCMEM");
        
    std::fstream fin; 
    std::fstream output_file("../data/output/testoutCMEM.dat"); // read output from birp
    std::fstream truth_file("../data/test/test_CMEM_static.dat"); // read output from a known good render
    std::fstream incorrect_file("../data/test/test_CMEM_static_incorrect.dat"); // just to be sure, try a known bad render

    std::vector<std::string> output_image;
    std::vector<std::string> truth_image;
    std::vector<std::string> incorrect_image;
    std::string output_row;
    std::string truth_row;
    std::string incorrect_row;

    while (std::getline(output_file, output_row, ',')) {
        output_image.push_back(output_row); // loop through and save all values
    }
    while (std::getline(truth_file, truth_row, ',')) {
        truth_image.push_back(truth_row); // loop through and save all values
    }
    while (std::getline(incorrect_file, incorrect_row, ',')) {
        incorrect_image.push_back(incorrect_row); // loop through and save all values
    }

    EXPECT_TRUE(std::equal( truth_image.begin(), truth_image.end(), output_image.begin() ));
    EXPECT_FALSE(std::equal( incorrect_image.begin(), incorrect_image.end(), output_image.begin() ));
}