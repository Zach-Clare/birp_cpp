#include <vector>
#include <cmath>
#include <memory>
#include <fstream>
#include <chrono>

#include "DataCube.h"
#include "Helper.h"
#include "Camera.h"

#include <EleFits/FitsFile.h>
#include <EleFits/SifFile.h>
#include <EleFits/MefFile.h>

Camera::Camera(DataCube cube, float pixel_size_deg, int plot_fov)
{
    ray_samples = 200;
    fov = plot_fov;

    dataCube = cube;
    image_dimension = std::round(plot_fov / pixel_size_deg);
    float image_data[image_dimension * image_dimension];
    // we woudl usually make the image array here, but we'll fill it dynamically

}

void Camera::SetPosition(float x, float y, float z)
{
    position.x = x;
    position.y = y;
    position.z = z;
}

void Camera::SetAim(float x, float y, float z)
{
    aim.x = x;
    aim.y = y;
    aim.z = z;
}

void Camera::Render()
{
    // this->Orient();
    // this->GetSky();
    this->GetRayStepDistances(44);
    this->Integrate();
}

void Camera::GetRayStepDistances(int max)
{
    for (int i = 0; i < ray_samples; i++) {
        ray_dist.push_back(max * ((i + 0.5f) / ray_samples));
    }

    ray_width = ray_dist[1] - ray_dist[0];

    return;
}

void Camera::Integrate() 
{
    float to_rad = M_PI / 180.f;
    float to_deg = 180.f / M_PI;

    // Let's calculate some angles. Firstly, a rotation about the z axis so that y points towards the aimpoint
    float x_diff = aim.x - position.x;
    float theta = (2 * M_PI) - std::atan2(x_diff, position.y);
    // float rotation_z = ((180 - (to_deg * theta)) * to_rad) + (20 * to_rad);
    float hypotenuse = std::sqrt(std::pow(x_diff, 2) + std::pow(position.y, 2));

    // This is wrong
    // float xh = std::sqrt(std::pow(position.y, 2) + std::pow(position.z, 2));
    // float angle = std::atan2(x_diff, xh);

    // Construct a right angle triangle the opposite side of the x-axis spanning from pos.x to aim.x
    // The hypotenuse of this triangle aligns on on the plane where x'z' intersects with xy. That hypotenuse is len.
    // float len = std::sqrt(std::pow(x_diff, 2) + std::pow(hypotenuse - position.y, 2));
    // float adjacent = std::sqrt(std::pow(hypotenuse, 2) + std::pow(position.z, 2));
    // float angle = std::atan2(len, adjacent);
    // Wrong again!

    // Third time's a charm
    float a = 90.f - (to_deg * std::atan2(position.y, x_diff));
    float far_y = x_diff * std::tan(a);
    float far_y_hypotenuse = std::sqrt(std::pow(far_y, 2) + std::pow(x_diff, 2));
    float total_y = std::abs(far_y) + std::abs(position.y);
    float back_hypotenuse = std::sqrt(std::pow(total_y, 2) + std::pow(position.z, 2));
    float angle = (2 * M_PI) - std::asin(far_y_hypotenuse / back_hypotenuse);
    
    float phi = std::atan2(hypotenuse, position.z);

    float rotation_z = 0.f * to_rad;
    std::vector<std::vector<float>> rz = {
        {std::cos(theta), - std::sin(theta), 0},
        {std::sin(theta), std::cos(theta), 0},
        {0, 0, 1}
    };
    std::vector<std::vector<float>> rx = {
        {1, 0, 0},
        {0, std::cos(phi), - std::sin(phi)},
        {0, std::sin(phi), std::cos(phi)}
    };
    // float angle = 5.f * to_rad;
    std::vector<std::vector<float>> ry = {
        {std::cos(angle), 0, std::sin(angle)},
        {0, 1, 0},
        {- std::sin(angle), 0, std::cos(angle)}
    };

    std::vector<std::vector<float>> identity = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    // Then we need to do matrix multiplication
    // https://en.wikipedia.org/wiki/Rotation_matrix#General_3D_rotations

    std::vector<std::vector<float>> rotation = Helper::MatrixMultiply(rx, rz);
    std::vector<std::vector<float>> rotation2 = Helper::MatrixMultiply(rotation, rx);
    std::vector<std::vector<float>> rotation_inverse = Helper::GetInverse(rotation);

    

    // serial method first, no parallelisation yet
    for (int i = 0; i < image_dimension; ++i) {
        for (int j = 0; j < image_dimension; ++j) {

            // Calculate pixel position
            float device_x = (j + 0.5f) / (image_dimension);
            float device_y = (i + 0.5f) / (image_dimension);
            float screen_x = ((2 * device_x) - 1) * std::tan(to_rad * (fov / 2));
            float screen_y = (1 - (2 * device_y)) * std::tan(to_rad * (fov / 2));
            float screen_vec[3] = {screen_x, screen_y, -1};

            float new_distance = Helper::VectorDistance(screen_vec);
            float screen_unit_vec[3] = {
                screen_x / new_distance,
                screen_y / new_distance,
                -1 / new_distance
            }; // That's our camera-based ray unit vector. We need to apply the rotation matrix
            
            std::vector<float> vector = Helper::ApplyRotation(rotation_inverse, screen_unit_vec);
            float* world_vector = &vector[0];
            float world_distance = Helper::VectorDistance(world_vector);
            std::vector<float> world_unit_vector = {
                world_vector[0] / world_distance,
                world_vector[1] / world_distance,
                world_vector[2] / world_distance
            }; // and now we have the unit vector in world space!

            int pxk_ok = 1;
            int pxk_yes = 0;

            for (int pxk = 0; pxk < pxn_dist; ++pxk) {
                if (pxk_ok == 1) {

                    // auto t_inner = std::chrono::high_resolution_clock::now();
                    
                    // vector to sample point
                    std::vector<float> sample_vector = Helper::MatrixScalarMultiply(world_unit_vector, ray_dist[pxk]);

                    // x coordinate of sample point
                    float x_coord = sample_vector[0] + position.x; // use nearest neighbour
                    float y_coord = sample_vector[1] + position.y;
                    float z_coord = sample_vector[2] + position.z;

                    // If this gives a segmentation fault, the datacube likely has a problem loading. TODO: Error catch
                    float x_ingress = std::abs(dataCube.coords_x[0] - x_coord);
                    int x_index = (x_ingress / dataCube.spacing[0]);

                    float y_ingress = std::abs(dataCube.coords_y[0] - y_coord);
                    int y_index = (y_ingress / dataCube.spacing[1]);

                    float z_ingress = std::abs(dataCube.coords_z[0] - z_coord);
                    int z_index = (z_ingress / dataCube.spacing[2]);

                    if ((0 > x_index || x_index >= dataCube.size.x) ||
                        (0 > y_index || y_index >= dataCube.size.y) ||
                        (0 > z_index || z_index >= dataCube.size.z)) {
                        pxk_ok = 0;
                        continue;
                    }

                    float sample = dataCube.slices.at(z_index).at(y_index).at(x_index);
                    sample = (sample * ray_width * (1000 * 100 * 6378.1) / (4 / M_PI)) / 10000;

                    // if (x_coord < (aim.x + 0.1f) && x_coord > (aim.x - 0.1f) &&
                    //     y_coord < 0.1f && y_coord > -0.1f &&
                    //     z_coord < 0.1f && z_coord > -0.1f) {
                    //         sample = 70;
                    // }

                    image[i][j] = image[i][j] + sample;
                    sample_vector.clear();

                    pxk_yes = 1;

                }
            }
        }
        // std::cout << std::to_string(i) << ", " << std::flush;
    }
}

int Camera::ToDat(std::string filename)
{
    std::ofstream outfile("../python/" + filename + ".dat");

    for (int i = 0; i < 144; ++i) {
        for (int j = 0; j < 144; ++j) {
            outfile << std::to_string(image[i][j]);
            outfile << ",";
        }
        outfile << "\n";
    }

    return 1;
}

int Camera::ToFITS(std::string filename)
{
    using namespace Euclid;

    std::vector<float> output_vector;
    for (int i = 0; i < image_dimension; i++) {
        for (int j = 0; j < image_dimension; j++) {
            output_vector.push_back(image[i][j]);
        }
    }

    float* outs = &output_vector[0];

    // // auto record = Fits::Record("rho", "1", "t", "emission from x-ray");
    long height = 144;
    long width = 144;
    auto raster = Fits::makeRaster(output_vector, height, width);
    Fits::SifFile f(filename + ".fits", Fits::FileMode::Create);
    // f.write("rho", "0", "t", "comment");
    // Fits::Record<std::string> record("rho", "0", "t", "comment");
    // f.write(record, raster);
    f.writeRaster(raster);
    
    // auto raster = Fits::makeRaster(outs, (long(144), long(144)));
    // Fits::MefFile f(filename, Fits::FileMode::Create);

    // const auto& image1 = f.append_image("IMAGE1", {}, raster);

    return 1;
}

Camera::~Camera()
{
}
