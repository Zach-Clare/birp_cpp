#include <cstdio>
#include <vector>
#include <string>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <unistd.h>

// #include "mainconfig.h"
#include "DataCube.h"
#include "Camera.h"
#include "Helper.h"

int main(int argc, char** argv)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();    

    int c;
    std::string input;
    std::string output;
    bool batch;

    // spacecraft position defaults
    float pos_x = 5.91813f;
    float pos_y = 6.74116f;
    float pos_z = 17.7067f;
    float aim = 7.81768f;

    while ((c = getopt(argc, argv, "i:o:x:y:z:a:")) != -1)
    {
        switch (c)
        {
            case 'i':
            {
                input = optarg;
                if (input.find(".dat") != std::string::npos) { // if arguements is a .dat file (a bit nieve)
                    batch = false;
                    std::cout << "non-batch";
                } else {
                    batch = true;
                    std::cout << "batch";
                }
                break;
            }
            case 'o':
            {
                // this should be a directory, but it might be worth checking
                output = optarg;
                std::cout << "output folder passed";
                break;
            }
            case 'x':
            {
                std::string fs(optarg);
                pos_x = std::stof(optarg);
                break;
            }
            case 'y':
            {
                std::string fs(optarg);
                pos_y = std::stof(optarg);
                break;
            }
            case 'z':
            {  
                std::string fs(optarg);
                pos_z = std::stof(optarg);
                break;
            }
            case 'a':
            {
                std::string fs(optarg);
                aim = std::stof(optarg);
                break;
            }
        }
    }


    DataCube cube;
    // cube.Load("/home/zc/code/birp/cpp/Batch/px_uni_0911.dat", 82, false);
    
    int plot_fov = 36;
    float pixel_size_deg = .25f;

    // Process in batch here
    if (batch) {
        // std::string batch_path = "/home/zc/code/birp/cpp/data/batch";
        for (const auto & entry : std::filesystem::directory_iterator(input)) {
            std::string path = entry.path().string();
            cube.Load(entry.path(), 82, false);

            Camera camera(cube, pixel_size_deg, plot_fov);
            camera.SetPosition(5.91813f, 6.74116f, 17.7067f);
            camera.SetAim(7.81768f, 0.f, 0.f);

            camera.Render();
            std::cout << "Rendered image.\n";
            std::string base_filename = path.substr(path.find_last_of("/\\") + 1);
            camera.ToFITS(output + base_filename);
            std::cout << "Exported.\n";
            // std::cout << entry.path() << std::endl;
        }
    } else {
        cube.Load(input, 82, false);

        Camera camera(cube, pixel_size_deg, plot_fov);
        camera.SetPosition(pos_x, pos_y, pos_z);
        camera.SetAim(aim, 0.f, 0.f);

        camera.Render();
        std::cout << "Rendered image.\n";
        std::string base_filename = input.substr(input.find_last_of("/\\") + 1);
        camera.ToFITS(output + base_filename);
        std::cout << "Exported.\n";
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    std::cout << duration.count() << std::endl;
}