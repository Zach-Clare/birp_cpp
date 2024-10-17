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

    // should the renderer use trilinear interpolation
    // default false. opt-in
    bool interpolate = false;

    while ((c = getopt(argc, argv, "i:o:x:y:z:a:t")) != -1)
    {
        switch (c)
        {
            case 'i': // input
            {
                input = optarg;
                if (input.find(".dat") != std::string::npos) { // if arguements is a .dat file (a bit nieve)
                    batch = false;
                    std::cout << "Using non-batch mode\n";
                } else {
                    batch = true;
                    std::cout << "Using batch mode\n";
                }
                break;
            }
            case 'o': // output
            {
                // this should be a directory, but it might be worth checking
                output = optarg;
                break;
            }
            case 'x': // x spacecraft pos
            {
                std::string fs(optarg);
                pos_x = std::stof(optarg);
                break;
            }
            case 'y': // y spacecraft pos
            {
                std::string fs(optarg);
                pos_y = std::stof(optarg);
                break;
            }
            case 'z': // z spacecraft pos
            {  
                std::string fs(optarg);
                pos_z = std::stof(optarg);
                break;
            }
            case 'a': // aim
            {
                std::string fs(optarg);
                aim = std::stof(optarg);
                break;
            }
            case 't': // trilinear - used as a flag
            {
                interpolate = true;
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
            std::cout << "Loading datacube..." << std::flush;
            cube.Load(entry.path(), 82, false);
            std::cout << "Datacube loaded.\nRendering..." << std::flush;

            Camera camera(cube, pixel_size_deg, plot_fov);
            camera.SetPosition(pos_x, pos_y, pos_z);
            camera.SetAim(aim, 0.f, 0.f);

            camera.Render(interpolate);
            std::cout << "Completed rendering.\nExporting..." << std::flush;
            std::string base_filename = path.substr(path.find_last_of("/\\") + 1);
            camera.ToFITS(output + base_filename);
            std::cout << "Exported." << std::flush;
            // std::cout << entry.path() << std::endl;
        }
    } else {
        std::cout << "Loading datacube...\t" << std::flush;
        cube.Load(input, 82, false);
        std::cout << "Datacube loaded.\nRendering...\t\t" << std::flush;

        Camera camera(cube, pixel_size_deg, plot_fov);
        camera.SetPosition(pos_x, pos_y, pos_z);
        camera.SetAim(aim, 0.f, 0.f);

        camera.Render(interpolate);
        std::cout << "Completed rendering.\nExporting...\t\t" << std::flush;
        std::string base_filename = input.substr(input.find_last_of("/\\") + 1);
        camera.ToFITS(output + base_filename);
        std::cout << "Exported.\n" << std::flush;
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    std::cout << "Entire operation took " << duration.count() << "ms" << std::endl;
}