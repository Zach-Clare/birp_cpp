#include <cstdio>
#include <vector>
#include <string>
#include <chrono>
#include <filesystem>

#include "mainconfig.h"
#include "Classes/DataCube.h"
#include "Classes/Camera.h"
#include "Helper.h"

int main(void)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();    

    DataCube cube;
    // cube.Load("/home/zc/code/birp/cpp/Batch/px_uni_0911.dat", 82, false);
    
    int plot_fov = 36;
    float pixel_size_deg = .25f;

    // Process in batch here
    std::string batch_path = "/home/zc/code/birp/cpp/Batch";
    for (const auto & entry : std::filesystem::directory_iterator(batch_path)) {
        std::string path = entry.path().string();
        cube.Load(entry.path(), 82, false);

        Camera camera(cube, pixel_size_deg, plot_fov);
        camera.SetPosition(5.91813f, 6.74116f, 17.7067f);
        camera.SetAim(7.81768f, 0.f, 0.f);

        camera.Render();
        std::cout << "Rendered image.\n";
        std::string base_filename = path.substr(path.find_last_of("/\\") + 1);
        camera.ToFITS(base_filename);
        std::cout << "Exported.\n";

        // std::cout << entry.path() << std::endl;
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    std::cout << duration.count();
}