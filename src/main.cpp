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
#include "Cmem.h"
#include "Space.h"

int main(int argc, char** argv)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();    

    int c;
    std::string input;
    std::string output;
    bool batch;
    bool use_cmem = false;

    // spacecraft position defaults
    float pos_x = 5.91813f;
    float pos_y = 6.74116f;
    float pos_z = 17.7067f;
    float aim = 7.81768f;
    // float aim = 4.f; // for debugging

    // // DEBUG
    // float pos_x = 5.91813f;
    // float pos_y = 20.f;
    // float pos_z = 0.f;
    // float aim = 5.f;

    // should the renderer use trilinear interpolation
    // default false. opt-in
    bool interpolate = false;

    while ((c = getopt(argc, argv, "i:o:x:y:z:a:t:c")) != -1)
    {
        switch (c)
        {
            case 'i': // input
            {
                input = optarg;
                if (input.find(".dat") != std::string::npos) { // if arguements is a .dat file (a bit nieve)
                    batch = false;
                    std::cout << "Using non-batch mode\n";
                } else if (Helper::ends_with(input, "/")) { // if arg ends with "/"
                    batch = true;
                    std::cout << "Using batch mode\n";
                } else { // no usable datacube given.
                    std::cout << "Please provide valid datacube or use -c to use CMEM." << std::endl;
                    exit(3);
                }
                break;
            }
            case 'c': // using CMEM
            {
                batch = false;
                use_cmem = true;
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

    // cube.Load("/home/zc/code/birp/cpp/Batch/px_uni_0911.dat", 82, false);
    
    // // Double fov - for debug and furthr context in renders
    // int plot_fov = 72;
    // float pixel_size_deg = .5f;

    // normal fov settings
    int plot_fov = 36;
    float pixel_size_deg = .25f;

    // Process in batch here
    if (batch) {
        // std::string batch_path = "/home/zc/code/birp/cpp/data/batch";
        for (const auto & entry : std::filesystem::directory_iterator(input)) {

            std::string path = entry.path().string();
            std::cout << "Loading datacube..." << std::flush;
            DataCube cube;
            cube.Load(entry.path(), 82, false);
            // space.Init(entry.path(), 82, false);
            std::cout << "Datacube loaded.\nRendering..." << std::flush;

            Camera camera(cube, pixel_size_deg, plot_fov);
            camera.SetPosition(pos_x, pos_y, pos_z);
            camera.SetAim(aim, 0.f, 0.f);

            camera.Render(interpolate);
            std::cout << "Completed rendering.\nExporting..." << std::flush;
            std::string base_filename = path.substr(path.find_last_of("/\\") + 1);
            camera.ToFITS(output + base_filename);
            std::cout << "Exported to " << output + base_filename << ".fits\n" << std::flush;
            // std::cout << entry.path() << std::endl;
        }
    } else { // non-batch

        std::string base_filename;

        Space *space; // Create a Space-shaped pointer
        if (!use_cmem) {
            std::string path = input;
            base_filename = input.substr(input.find_last_of("/\\") + 1);
            std::cout << "Loading datacube...\t" << std::flush;
            DataCube* cube = new DataCube(); // Create the datacube object and get the pointer
            space = cube; // Save the pointer to DataCube in our Space-shaped hole
            cube->Init(path, 82, false); // Use as cube, not space (due to Object Slicing)
            cube->SetTrilinear(interpolate); // save command-line option to object
            std::cout << "Datacube loaded.\nRendering...\t\t" << std::flush;
        } else {
            // this is where CMEM is used - no need to load data
            base_filename = "cmembirp";
            std::cout << "Initialising CMEM...\t" << std::flush;
            CMEM* cmem = new CMEM();
            space = cmem;
            cmem->Init();
            std::cout << "Initialised.\nRendering...\t\t" << std::flush;
        }

        Camera camera(*space, pixel_size_deg, plot_fov); // Init camera with pointer to space rather than copying the whole object over
        camera.SetPosition(pos_x, pos_y, pos_z);
        camera.SetAim(aim, 0.f, 8.f);

        camera.Render(interpolate);
        delete space;
        std::cout << "Completed rendering.\nExporting...\t\t" << std::flush;
        camera.ToFITS(output + base_filename);
        std::cout << "Exported to " << output << base_filename << ".fits\n" << std::flush;
    }

    auto end = high_resolution_clock::now();
    auto duration_ms = duration_cast<milliseconds>(end - start);
    auto duration_s = duration_cast<seconds>(end - start);
    std::cout << "Entire operation took " << duration_ms.count() << "ms or " << duration_s.count() << "s" << std::endl;
}