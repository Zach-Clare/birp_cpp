#include <iostream>
#include <vector>
#include <string>
// #include <EleFits/MefFile.h>

#include "Classes/DataCube.h"
#include "Classes/Camera.h"
#include "Helper.h"

int main()
{
    // Here, we'll use as a substitute for main.py

    DataCube cube;
    cube.Load("Data/px_uni_1000.dat", 82, false);

    int plot_fov = 36;
    float pixel_size_deg = .25f;

    Camera camera(cube, pixel_size_deg, plot_fov);
    camera.SetPosition(5.91813f, 6.74116f, 17.7067f);
    camera.SetAim(7.81768f, 0.f, 0.f);

    // camera.SetPosition(7.70111f, 7.86787f, 9.91106f);
    // camera.SetAim(6.06733f, 0.f, 0.f);

    camera.Render();
    camera.ToFITS();
}