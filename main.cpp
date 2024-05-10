#include <iostream>
#include <vector>
#include <string>
#include "Classes/DataCube.h"
#include "Classes/Camera.h"

int main()
{
    // Here, we'll use as a substitute for main.py

    DataCube cube;
    cube.Load("Data/px_uni_1000.dat", 82);

    int plot_fov = 36;
    float pixel_size_deg = 0.25f;

    Camera camera(cube, pixel_size_deg, plot_fov);

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