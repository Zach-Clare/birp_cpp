#include <iostream>
#include <vector>
#include <string>
#include "Classes/DataCube.cpp"

int main()
{
    // Here, we'll use as a substitute for main.py

    DataCube cube;
    cube.load("Data/px_uni_1000.dat", 82);
}

void amr_ov_cube2image(
    float orb[3],
    float aim[3], 
    float imgWidthDeg = 50.0,
    float pixSizeDeg = 1.0,
    std::string pxEmiss = "path/to/default/data/cube",
    int pxNDist = 90,
    float prMax = 45.0,
    std::string outroot = "pxout",
    std::string pxinGSM = "N"
    /// more to come, we'll do the rest as they appear
    )
{
    // read pxEmiss file
}