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

    // // spacecraft position defaults
    float pos_x = 5.91813f;
    float pos_y = 6.74116f;
    float pos_z = 17.7067f;
    float aim = 7.81768f;
    // float aim = 4.f; // for debugging

    // float pos_x = -1.07f;
    // float pos_y = 0.26f;
    // float pos_z = 12.42f;
    // float aim = 6.05f;

    // // DEBUG
    // float pos_x = 5.91813f;
    // float pos_y = 20.f;
    // float pos_z = 0.f;
    // float aim = 5.f;

    // should the renderer use trilinear interpolation
    // default false. opt-in
    bool interpolate = false;
    
    // camera properties defaults
    float pixel_size_deg = .25f;
    // int plot_fov = 108; // Double fov - for debug and furthr context in renders
    // float pixel_size_deg = .75f; // double
    float plot_fov_h = 36.f; // normal fov settings
    float plot_fov_w = 36.f; // normal

    // Init CMEM params 
    std::vector<int> v; // int(3)
    std::vector<int> b; // int(3)
    float dipole = NULL; // dipole tilt
    std::vector<float> p = {1.f, 1.f, 3.f, 4.f}; // float (3)
    float B = NULL; // captial B param
    float alpha = NULL; // They should name these differently
    float beta = NULL;  // because that's a bit confusing
    float bs = NULL; // stands for Bow Shock, location of bow shock
    float A1 = NULL; // parameter values, unsure what these do
    float A2 = NULL;
    float ay_bs = NULL; // flaring parameters for bowshock
    float az_bs = NULL;

    while ((c = getopt(argc, argv, "i:o:x:y:z:a:tcs:h:w:v:b:d:p:q:r:f:e:g:u:j:k:")) != -1) //JKLMN remain. Time for Boost.program_options?
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
                    std::cout << "Please provide valid .dat datacube or use -c to use CMEM." << std::endl;
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
            case 's': // pixel size in degrees
            {
                std::string fs(optarg);
                pixel_size_deg = std::stof(optarg);
                break;
            }
            case 'h': // height FOV of image
            {
                std::string fs(optarg);
                plot_fov_h = std::stof(optarg);
                break;
            }
            case 'w': // width FOV of image
            {
                std::string fs(optarg);
                plot_fov_w = std::stof(optarg);
                break;
            }
            case 'v': // v solar wind param for CMEM ///////////////////////////////////////////////////// unused?
            {
                std::string fs(optarg);
                v = Helper::explode_int(fs, ','); // process, split by comma
                break;
            }
            case 'b': // b solar wind param for CMEM ////////////////////////////////////////////////////// unused?
            {
                std::string fs(optarg);
                b = Helper::explode_int(fs, ','); // process, split by comma
                break;
            }
            case 'd': // dipole tilt for CMEM ///////////////////////////////////////////////////////////// unused?
            {
                std::string fs(optarg);
                dipole = std::stof(fs);
                break;
            }
            case 'p': // (p0, p1, p2, p3) solar wind param for CMEM // scaling factor for bs distance, mp flaring, and mp indentation, and unknown
            {
                std::string fs(optarg);
                p = Helper::explode_float(fs, ','); // process, split by comma
                break;
            }
            case 'q': // B running out of letters //  CMEM-specific maps to B param "big B"
            {
                std::string fs(optarg);
                B = std::stof(fs);
                break;
            } 
            case 'r': // alpha // CMEM-specific
            {
                std::string fs(optarg);
                alpha = std::stof(fs);
                break;
            }
            case 'f': // beta // CMEM-specific
            {  
                std::string fs(optarg);
                beta = std::stof(fs);
                break;
            }
            case 'e': // bs bowshock distance // CMEM-specific
            {
                std::string fs(optarg);
                bs = std::stof(fs);
                break;
            }
            case 'g': // A1 CMEM-specific
            {
                std::string fs(optarg);
                A1 = std::stof(fs);
                break;
            }
            case 'u': // A2 CMEM-specific
            {
                std::string fs(optarg);
                A2 = std::stof(fs);
                break;
            }
            case 'j': // bs flaring parameter in y CMEM-specific
            {
                std::string fs(optarg);
                ay_bs = std::stof(fs);
                break;
            }
            case 'k': // bs flaring parameter in z CMEM-specific
            {
                std::string fs(optarg);
                az_bs = std::stof(fs);
                break;
            }
        }
    }

    // cube.Load("/home/zc/code/birp/cpp/Batch/px_uni_0911.dat", 82, false);

    // Process in batch here
    if (batch) {
        // std::string batch_path = "/home/zc/code/birp/cpp/data/batch";
        for (const auto & entry : std::filesystem::directory_iterator(input)) {

            std::string path = entry.path().string();
            std::cout << "Loading datacube..." << std::flush;
            DataCube cube;
            cube.Load(entry.path(), 82, false);
            cube.SetTrilinear(interpolate);
            // space.Init(entry.path(), 82, false);
            std::cout << "Datacube loaded.\nRendering..." << std::flush;

            Camera camera(cube, pixel_size_deg, plot_fov_h, plot_fov_w);
            camera.SetPosition(pos_x, pos_y, pos_z);
            camera.SetAim(aim, 0.f, 0.f);

            camera.Render();
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
            if (dipole != NULL ||
                B != NULL ||
                alpha != NULL ||
                beta != NULL || 
                !v.empty() ||
                !b.empty() ||
                !p.empty()
            ) {
                cmem->Init(true, v, b, dipole, p[0], p[1], p[2], p[3], B, alpha, beta, bs, A1, A2, ay_bs, az_bs);
            } else {
                cmem->Init();
            }
            std::cout << "Initialised.\nRendering...\t\t" << std::flush;
        }

        Camera camera(*space, pixel_size_deg, plot_fov_h, plot_fov_w); // Init camera with pointer to space rather than copying the whole object over
        camera.SetPosition(pos_x, pos_y, pos_z);
        camera.SetAim(aim, 0.f, 0.f);

        camera.Render();
        // delete space;
        std::cout << "Completed rendering.\nExporting...\t\t" << std::flush;
        camera.ToFITS(output + base_filename);
        // camera.ToDat(output + base_filename);
        std::cout << "Exported to " << output << base_filename << ".fits\n" << std::flush;
    }

    auto end = high_resolution_clock::now();
    auto duration_ms = duration_cast<milliseconds>(end - start);
    auto duration_s = duration_cast<seconds>(end - start);
    std::cout << "Entire operation took " << duration_ms.count() << "ms or " << duration_s.count() << "s" << std::endl;
}