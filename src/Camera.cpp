/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


#include <vector>
#include <cmath>
#include <memory>
#include <fstream>
#include <chrono>
#include <filesystem>

#include "DataCube.h"
#include "Cmem.h"
#include "Space.h"
#include "Helper.h"
#include "Camera.h"

#include <EleFits/FitsFile.h>
#include <EleFits/SifFile.h>
#include <EleFits/MefFile.h>
#include <EleFits/Header.h>
#include <EleFitsData/Raster.h>

Camera::Camera(Space &cube, float pixel_size_degrees, float plot_fov_h, float plot_fov_w)
{
    ray_samples = 200; // now many samples to take along each ray
    fov_x = plot_fov_w;
    fov_y = plot_fov_h;
    pixel_size_rad = pixel_size_deg * (M_PI / 180.f);
    pixel_size_deg = pixel_size_degrees;

    dataCube = &cube;
    image_dimension_y = std::round(plot_fov_h / pixel_size_deg);
    image_dimension_x = std::round(plot_fov_w / pixel_size_deg);
    // image_dimension_x = image_x;
    // image_dimension_y = image_y;
    // float image_data[image_dimension_x * image_dimension_y];
    // we woudl usually make the image array here, but we'll fill it dynamically

    std::vector<std::vector<float>> image_whatevs(image_dimension_y, std::vector<float>(image_dimension_x, 0.f));
    image = image_whatevs;
}

void Camera::SetPosition(float& x, float& y, float& z)
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
    this->GenerateRayDistWidth(44);
    this->Integrate();
}

void Camera::GenerateRayDistWidth(int max)
{
    float dist;
    for (int i = 0; i < ray_samples; i++) {
        dist = max * ((i + 0.5f) / ray_samples);
        ray_dist.push_back(dist);
        ray_widths.push_back((sin(pixel_size_rad / 2) * dist) * 2);
    }

    ray_width = ray_dist[1] - ray_dist[0];

    return;
}

void Camera::Integrate() 
{
    float to_rad = M_PI / 180.f;
    // float to_deg = 180.f / M_PI; // uncomment if using

    float angle_sun = Camera::Orient();
    float angle_sun_corrected = 180.f - angle_sun; /// CHANGED 360 TO 180

    // Let's calculate some angles. Firstly, a rotation about the z axis so that y points towards the aimpoint
    float x_diff = aim.x - position.x;
    float theta = (2 * M_PI) - std::atan2(x_diff, position.y);
    // float rotation_z = ((180 - (to_deg * theta)) * to_rad) + (20 * to_rad);
    float hypotenuse = std::sqrt(std::pow(x_diff, 2) + std::pow(position.y, 2));

    // Then calculate the rotation about the x axis
    float phi = std::atan2(hypotenuse, position.z);

    // float rotation_z = 0.f * to_rad; // would be used in rz matrix
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
    std::vector<std::vector<float>> rz2 = {
        {std::cos(angle_sun_corrected * to_rad), - std::sin(angle_sun_corrected * to_rad), 0},
        {std::sin(angle_sun_corrected * to_rad), std::cos(angle_sun_corrected * to_rad), 0},
        {0, 0, 1}
    };

    // float angle = 5.f * to_rad;
    // std::vector<std::vector<float>> ry = {
    //     {std::cos(angle), 0, std::sin(angle)},
    //     {0, 1, 0},
    //     {- std::sin(angle), 0, std::cos(angle)}
    // };

    std::vector<std::vector<float>> identity = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    // Then we need to do matrix multiplication
    // https://en.wikipedia.org/wiki/Rotation_matrix#General_3D_rotations
    
    std::vector<std::vector<float>> rotation = Helper::MatrixMultiply(rz2, rx);
    std::vector<std::vector<float>> rotation2 = Helper::MatrixMultiply(rotation, rz);
    std::vector<std::vector<float>> rotation_inverse = Helper::GetInverse(rotation2);

    float aspect_ratio = (float)image_dimension_x / (float)image_dimension_y;

    // serial method first, no parallelisation yet
    for (int j = 0; j < image_dimension_y; ++j) {
        for (int i = 0; i < image_dimension_x; ++i) {

            // i and j are raster space (total number of pixels)
            // device_x and _y are NDC space (normalised from 0 - 1)
            // screen_x and _y are screen space (normalised from -1 - 1)

            // Calculate pixel position
            float device_x = (i + 0.5f) / (image_dimension_x);
            float device_y = (j + 0.5f) / (image_dimension_y);
            float screen_x = ((2 * device_x) - 1) * std::tan(to_rad * (fov_x / 2)) * aspect_ratio;
            float screen_y = (1 - (2 * device_y)) * std::tan(to_rad * (fov_y / 2));
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
            // int pxk_yes = 0;

            for (int pxk = 0; pxk < ray_samples; ++pxk) {
                if (pxk_ok == 1) {

                    // auto t_inner = std::chrono::high_resolution_clock::now();
                    
                    // // Reverse ray direction, start at the back
                    // This is to debug the rendering process
                    // float debug_dist =  44 - ray_dist[pxk];

                    // vector to sample point
                    std::vector<float> sample_vector = Helper::MatrixScalarMultiply(world_unit_vector, ray_dist[pxk]);
                    // std::vector<float> sample_vector = Helper::MatrixScalarMultiply(world_unit_vector, debug_dist); // if wanting to reverse

                    // x coordinate of sample point
                    float x_coord = sample_vector[0] + position.x; // use nearest neighbour
                    float y_coord = sample_vector[1] + position.y;
                    float z_coord = sample_vector[2] + position.z;

                    float sample = dataCube->GetSample(x_coord, y_coord, z_coord);

                    // float width_factor;
                    // if (pxk == 0) {
                    //     width_factor = ray_widths[pxk] * ray_widths[pxk] * ray_dist[pxk];
                    // } else {
                    //     width_factor = ray_widths[pxk] * ray_widths[pxk] * (ray_dist[pxk] - ray_dist[pxk - 1]);
                    // }

                    // sample = (sample * width_factor * (1000 * 100 * 6378.1) / (4 / M_PI)) / 10000;
	                sample = (sample * ray_width * (1000 * 100 * 6378.1) / (4 / M_PI)) / 10000;
                    // sample = (sample * ray_width);
                    // sample = sample * 3;

                    // if sample is negative, location is out of bounds
                    // ignore remainder of this ray
                    if (sample < 0) {
                        pxk_ok = 0;
                        continue;
                    }

                    // //// DEBUG ONLY - this shows the x line in the final render
                    // if (y_coord < 0.1f && y_coord > -0.1f &&
                    //     z_coord < 0.1f && z_coord > -0.1f) {
                    //         sample = 70;
                    // }

                    // if (i == 0 && j == 0) {
                    //     sample = 0.1f;
                    // }

                    // image[i][j] = image[i][j] + sample;
                    image[j][i] = image[j][i] + sample; // one of these needs to be commented out?
                    sample_vector.clear();

                    // pxk_yes = 1;

                }
            }
        }
        // std::cout << std::to_string(i) << ", " << std::flush;
    }
}

int Camera::ToDat(std::string filename)
{
    std::ofstream outfile(filename + ".dat");

    for (int j = 0; j < image_dimension_y; ++j) {
        for (int i = 0; i < image_dimension_x; ++i) {
            outfile << std::to_string(image[j][i]);
            outfile << ",";
        }
        outfile << "\n";
    }

    outfile.close();

    return 1;
}

int Camera::ToFITS(std::string filename)
{
    using namespace Euclid;

    std::vector<float> output_vector;
    for (int j = 0; j < image_dimension_y; j++) {
        for (int i = 0; i < image_dimension_x; i++) {
            output_vector.push_back(image[j][i]);
        }
    }

    auto raster = Fits::makeRaster(std::move(output_vector), image_dimension_x, image_dimension_y);
    Fits::Record<std::string> record("rho", "0", "t", "comment");
    Fits::Record<float> crval1 {"CRVAL1", -(this->fov_x / 2), "deg", ""};
    Fits::Record<float> cdelt1 {"CDELT1", pixel_size_deg, "deg", ""};
    Fits::Record<float> crval2 {"CRVAL2", -(this->fov_y / 2), "deg", ""};
    Fits::Record<float> cdelt2 {"CDELT2", pixel_size_deg, "deg", ""};

    
    try {
        Fits::SifFile f(filename + ".fits", Fits::FileMode::Create);
        // f.write(crval1, raster);
        // f.write(cdelt1, raster);
        // f.write(crval2, raster);
        // f.write(cdelt2, raster);
        // Fits::Header()
        f.header().writeSeq(crval1, cdelt1, crval2, cdelt2);
        f.write(record, raster);

    } catch (Euclid::Cfitsio::CfitsioError) {
        // Error here is either an error with EleFits or the file exists already
        // So let's attempt to delete an existing file and try again
        // If the file does not exist already, there's a problem with EleFits, good luck

        std::cout << "\x1B[31mOverwriting\033[0m - " << std::flush; // This may not work with windows
        std::filesystem::remove(filename);
        Fits::SifFile f(filename + ".fits", Fits::FileMode::Overwrite);
        // f.write(crval1, raster);
        // f.write(cdelt1, raster);
        // f.write(crval2, raster);
        // f.write(cdelt2, raster);
        
        f.header().writeSeq(crval1, cdelt1, crval2, cdelt2);
        f.write(record, raster);

    }

    return 1;
}

float Camera::Orient()
{
    // This is a port of the IDL version of orient() so I'm
    // not entirely sure what it's doing, but I think it finds
    // out how to rotate the camera so that the x-line is 
    // perpendicular to the sides of the image.

    // I was going to refactor this, but it takes less than 1ms.
    // It's just maths with no loops so it's quick

    float to_rad = M_PI / 180.f;
    float to_deg = 180.f / M_PI;

    // get unit vector
    std::vector<float> distance_vec {
        aim.x - position.x,
        aim.y - position.y,
        aim.z - position.z
    };

    float distance = Helper::VectorDistance(&distance_vec[0]);
    std::vector<float> unit_vec = {
        distance_vec[0] / distance,
        distance_vec[1] / distance,
        distance_vec[2] / distance
    };

    // now we calculate the north vector
    float phi = std::acos(unit_vec[2]);
    float denominator = std::sin(phi);

    float theta1;
    if (std::abs(denominator) < std::abs(unit_vec[0])) {
        if (denominator * unit_vec[0] < 0) {
            theta1 = std::acos(-1) * to_deg;
        } else {
            theta1 = std::acos(1) * to_deg;
        }
    } else {
        theta1 = std::asin(unit_vec[0] / denominator);
    }

    float theta2;
    if (std::abs(denominator) < std::abs(unit_vec[1])) {
        if (denominator * unit_vec[1] < 0) {
            theta2 = std::acos(-1) * to_deg;
        } else {
            theta2 = std::acos(1) * to_deg;
        }
    } else {
        theta2 = std::asin(unit_vec[1] / denominator);
    }

    std::vector<float> north = {
        -1.f * std::cos(theta2) * std::cos(phi),
        -1.f * std::sin(theta2) * std::cos(phi),
        std::sin(phi)
    };

    if (unit_vec[0] >= 0 && unit_vec[1] >= 0 && unit_vec[2] >= 0) {
        north[0] = -1.f * std::abs(north[0]);
        north[1] = -1.f * std::abs(north[1]);
    } else if (unit_vec[0] >= 0 && unit_vec[1] >= 0 && unit_vec[2] < 0) {
        north[0] = 1.f * std::abs(north[0]);
        north[1] = 1.f * std::abs(north[1]);
    } else if (unit_vec[0] >= 0 && unit_vec[1] < 0 && unit_vec[2] >= 0) {
        north[0] = -1.f * std::abs(north[0]);
        north[1] = 1.f * std::abs(north[1]);
    } else if (unit_vec[0] >= 0 && unit_vec[1] < 0 && unit_vec[2] < 0) {
        north[0] = 1.f * std::abs(north[0]);
        north[1] = -1.f * std::abs(north[1]);
    } else if (unit_vec[0] < 0 && unit_vec[1] >= 0 && unit_vec[2] >= 0) {
        north[0] = 1.f * std::abs(north[0]);
        north[1] = -1.f * std::abs(north[1]);
    } else if (unit_vec[0] < 0 && unit_vec[1] >= 0 && unit_vec[2] < 0) {
        north[0] = -1.f * std::abs(north[0]);
        north[1] = 1.f * std::abs(north[1]);
    } else if (unit_vec[0] < 0 && unit_vec[1] < 0 && unit_vec[2] >= 0) {
        north[0] = 1.f * std::abs(north[0]);
        north[1] = 1.f * std::abs(north[1]);
    } else if (unit_vec[0] < 0 && unit_vec[1] < 0 && unit_vec[2] < 0) {
        north[0] = -1.f * std::abs(north[0]);
        north[1] = -1.f * std::abs(north[1]);
    }

    north[2] = std::abs(north[2]);

    // calculate the right vector
    std::vector<float> right = {
        (unit_vec[1] * north[2]) - (unit_vec[2] * north[1]),
        (unit_vec[2] * north[0]) - (unit_vec[0] * north[2]),
        (unit_vec[0] * north[1]) - (unit_vec[1] * north[0])
    };

    float plac = (-1.f * distance_vec[0] * position.x) + (-1.f * distance_vec[1] * position.y) + (-1.f * distance_vec[2] * position.z);

    std::vector<float> sun = {23455.f, 0.f, 0.f};
    std::vector<float> p2pout = PointToPlane(sun, plac, 1, distance_vec, unit_vec, north, right);
    float angle_sun = 270 - (to_deg * std::atan2(0 - p2pout[0], 0 - p2pout[1]));

    if (angle_sun > 360) {
        angle_sun = angle_sun - 360;
    } else if (angle_sun < 0) {
        angle_sun = angle_sun + 360;
    }

    return angle_sun;

}

std::vector<float> Camera::PointToPlane(std::vector<float> obj, float plac, float prad, std::vector<float> distance_vec, std::vector<float> unit_vec, std::vector<float> north, std::vector<float> right)
{
    float to_deg = 180 / M_PI;

    float snumer = (distance_vec[0] * obj[0]) + (distance_vec[1] * obj[1]) + (distance_vec[2] * obj[2]) + plac;
    float sdenom = (distance_vec[0] * unit_vec[0]) + (distance_vec[1] * unit_vec[1]) + (distance_vec[2] * unit_vec[2]);
    float scalar = snumer / sdenom;

    float q[3] = {
        obj[0] - (scalar * unit_vec[0]),
        obj[1] - (scalar * unit_vec[1]),
        obj[2] - (scalar * unit_vec[2])
    };

    // vector of object from the camera
    float object_vector[3] = {
        q[0] - position.x,
        q[1] - position.y,
        q[2] - position.z
    };

    float dp = std::sqrt((object_vector[0] * object_vector[0]) + (object_vector[1] * object_vector[1]) + (object_vector[2] * object_vector[2]));
    float object_unit_vector[3] = {
        object_vector[0] / dp,
        object_vector[1] / dp,
        object_vector[2] / dp
    };
    dp = std::sqrt((object_unit_vector[0] * object_unit_vector[0]) + (object_unit_vector[1] * object_unit_vector[1]) + (object_unit_vector[2] * object_unit_vector[2]));

    // angle in plane, angle of object from the north vector
    float t1 = (object_unit_vector[0] * north[0]) + (object_unit_vector[1] * north[1]) + (object_unit_vector[2] * north[2]);
    float t2 = dp * std::sqrt((north[0] * north[0]) + (north[1] * north[1]) + (north[2] * north[2]));
    float t3 = t1 / t2;
    if (t3 > 1) {
        t3 = 1.f;
    } else if (t3 < -1) {
        t3 = -1.f;
    }
    float apon = to_deg * std::acos(t3); // not sure what apon means (angle in plane of north?)

    // now get the angle of object from right vector
    float t4 = (object_unit_vector[0] * right[0]) + (object_unit_vector[1] * right[1]) + (object_unit_vector[2] * right[2]);
    float t5 = dp * std::sqrt((right[0] * right[0]) + (right[1] * right[1]) + (right[2] * right[2]));
    float t6 = t4/t5;
    if (t6 > 1) {
        t6 = 1.f;
    } else if (t6 < -1) {
        t6 = -1.f;
    }
    float apor = to_deg * std::acos(t6);
    float soy = 1.f;
    float sox = 1.f;
    if (apor > 90) {
        sox = -1.f;
    }

    // Vector from observer to object
    float observer_object_vector[3] = {
        obj[0] - position.x,
        obj[1] - position.y,
        obj[2] - position.z
    };
    float observer_object_distance = std::sqrt((observer_object_vector[0] * observer_object_vector[0]) + (observer_object_vector[1] * observer_object_vector[1]) + (observer_object_vector[2] * observer_object_vector[2]));
    // turn vector into unit vector
    observer_object_vector[0] = observer_object_vector[0] / observer_object_distance;
    observer_object_vector[1] = observer_object_vector[1] / observer_object_distance;
    observer_object_vector[2] = observer_object_vector[2] / observer_object_distance;

    // angular radius object (radius of the object)
    float aro = to_deg * std::atan(prad / observer_object_distance);

    // angle of observer object aim
    float t7 = (observer_object_vector[0] * unit_vec[0]) + (observer_object_vector[1] * unit_vec[1]) + (observer_object_vector[2] * unit_vec[2]);
    float t8 = std::sqrt((observer_object_vector[0] * observer_object_vector[0]) + (observer_object_vector[1] * observer_object_vector[1]) + (observer_object_vector[2] * observer_object_vector[2]));
    float t9 = t7 / t8;
    if (t9 > 1) {
        t9 = 1.f;
    }
    float observer_object_angle = to_deg * std::acos(t9);

    float xo = 0.f;
    float yo = 0.f;
    if (observer_object_angle > 0.000001f) {
        xo = sox * observer_object_angle * std::sin(apon / to_deg);
        yo = soy * observer_object_angle * std::cos(apon / to_deg);
    }

    std::vector<float> result = {xo, yo, aro, observer_object_angle};
    return result;
}

Camera::~Camera()
{
}
