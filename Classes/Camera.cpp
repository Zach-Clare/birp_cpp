#include <vector>
#include <cmath>
#include <memory>
#include <fstream>
#include <chrono>

#include "DataCube.h"
#include "../Helper.h"

#include "Camera.h"

Camera::Camera(DataCube cube, float pixel_size_deg, int plot_fov)
{
    ray_samples = 200;
    fov = plot_fov;

    dataCube = cube;
    image_dimension = std::round(plot_fov / pixel_size_deg);
    float image_data[image_dimension * image_dimension];
    // we woudl usually make the image array here, but we'll fill it dynamically

    this->BuildLatLon(pixel_size_deg, plot_fov);

    std::cout << "Lat and Lon arrays complete.\n";

    int n_sky_30 = 10;
    skysep = 30 / n_sky_30;
    ns_lat = std::floor(180 / skysep) + 1;
    ns_lon = std::floor(360 / skysep);
    nsky = ns_lat * ns_lon;

    this->GenerateSky(pixel_size_deg);

    std::cout << "Sky generation complete.\n";

    gei_to_gse = Helper::Gei2GseTransforms();

    std::cout << "Transforms complete.\n";

}

// Build the lat and lon arrays. Yep, with identical data.
void Camera::BuildLatLon(float pixel_size_deg, int plot_fov)
{
    float factor = (0.5f * pixel_size_deg) - (0.5f * plot_fov);
    for (int i = 0; i < image_dimension; i++) {
        lat.push_back(factor + (pixel_size_deg * i));
        lon.push_back(factor + (pixel_size_deg * i));
    }
}

void Camera::GenerateSky(float pixel_size_deg)
{
    int c_lat = 1;
    int c_lon = 0;

    for (int i = 0; i < nsky; i++) {
        if ((i + 1) > (c_lat * ns_lon)) {
            c_lat++;
        }

        float value = 90.0f - ((c_lat - 1) * skysep);
        if (value == 90.0f) {
            value = std::max(90.0f - (pixel_size_deg / 2), 89.5f);
        } else if (value == -90.0) {
            value = std::min(-90.0f - (pixel_size_deg / 2), -89.5f);
        }
        sky2.push_back(value);
        sky1.push_back(c_lon * skysep);

        c_lon++;

        if (c_lon == ns_lon) {
            c_lon = 0;
        }


        // Fill out the xsky and ysky vectors
        xsky.push_back(0.0f);
        ysky.push_back(0.0f);
    }

    // not sure why we do this. It was in the original IDL program
    xsky[nsky - 1] = -999.9f;
    ysky[nsky - 1] = -999.9f;
    
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
    this->Orient();
    this->GetSky();
    this->GetRayStepDistances(33);
    this->Integrate();
}

std::vector<float> Camera::PointToPlane(std::vector<float> obj, float plac, float prad)
{
    float conv = 180 / M_PI;

    float snumer = (dist[0] * obj[0]) + (dist[1] * obj[1]) + (dist[2] * obj[2]) + plac;
    float sdenom = (dist[0] * aimpoint_unit_vector[0]) + (dist[1] * aimpoint_unit_vector[1]) + (dist[2] * aimpoint_unit_vector[2]);
    float scalar = snumer / sdenom;

    float q[3] = {
        obj[0] - (scalar * aimpoint_unit_vector[0]),
        obj[1] - (scalar * aimpoint_unit_vector[1]),
        obj[2] - (scalar * aimpoint_unit_vector[2])
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
    float t1 = (object_unit_vector[0] * north_vector[0]) + (object_unit_vector[1] * north_vector[1]) + (object_unit_vector[2] * north_vector[2]);
    float t2 = dp * std::sqrt((north_vector[0] * north_vector[0]) + (north_vector[1] * north_vector[1]) + (north_vector[2] * north_vector[2]));
    float t3 = t1 / t2;
    if (t3 > 1) {
        t3 = 1.f;
    } else if (t3 < -1) {
        t3 = -1.f;
    }
    float apon = conv * std::acos(t3); // not sure what apon means (angle in plane of north?)

    // now get the angle of object from right vector
    t1 = (object_unit_vector[0] * right_vector[0]) + (object_unit_vector[1] * right_vector[1]) + (object_unit_vector[2] * right_vector[2]);
    t2 = dp * std::sqrt((right_vector[0] * right_vector[0]) + (right_vector[1] * right_vector[1]) + (right_vector[2] * right_vector[2]));
    t3 = t1/t2;
    if (t3 > 1) {
        t3 = 1.f;
    } else if (t3 < -1) {
        t3 = -1.f;
    }
    float apor = conv * std::acos(t3);
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
    float aro = conv * std::atan(prad / observer_object_distance);

    // angle of observer object aim
    t1 = (observer_object_vector[0] * aimpoint_unit_vector[0]) + (observer_object_vector[1] * aimpoint_unit_vector[1]) + (observer_object_vector[2] * aimpoint_unit_vector[2]);
    t2 = std::sqrt((observer_object_vector[0] * observer_object_vector[0]) + (observer_object_vector[1] * observer_object_vector[1]) + (observer_object_vector[2] * observer_object_vector[2]));
    t3 = t1 / t2;
    if (t3 > 1) {
        t3 = 1.f;
    }
    float observer_object_angle = conv * std::acos(t3);

    float xo = 0.f;
    float yo = 0.f;
    if (observer_object_angle > 0.000001f) {
        xo = sox * observer_object_angle * std::sin(apon / conv);
        yo = soy * observer_object_angle * std::cos(apon / conv);
    }

    std::vector<float> result = {xo, yo, aro, observer_object_angle};
    return result;
}

void Camera::Orient()
{
    float conv = 180 / M_PI;

    dist[0] = aim.x - position.x;
    dist[1] = aim.y - position.y;
    dist[2] = aim.z - position.z;

    float distance = std::sqrt((std::pow(dist[0], 2)) + (std::pow(dist[1], 2)) + (std::pow(dist[2], 2)));

    // aimpoint vector as a unit vector
    aimpoint_unit_vector[0] = dist[0] / distance;
    aimpoint_unit_vector[1] = dist[1] / distance;
    aimpoint_unit_vector[2] = dist[2] / distance;
    
    // calculate the north vector
    float phi = conv * std::acos(aimpoint_unit_vector[2]);
    float thdenom = std::sin(phi / conv);

    // calculate based around unit vectors x value
    float theta1;
    if (std::abs(thdenom) < std::abs(aimpoint_unit_vector[0])) {
        if ((thdenom * aimpoint_unit_vector[0]) < 0) {
            theta1 = conv * std::acos(-1);
        } else if ((thdenom * aimpoint_unit_vector[0]) > 0) {
            theta1 = conv * std::acos(1);
        }
    } else {
        theta1 = conv * std::acos(aimpoint_unit_vector[0] / thdenom);
    }

    // calculate based around unit vectors y value
    float theta2;
    if (std::abs(thdenom) < std::abs(aimpoint_unit_vector[1])) {
        if ((thdenom * aimpoint_unit_vector[1]) < 0) {
            theta2 = conv * std::acos(-1);
        } else if ((thdenom * aimpoint_unit_vector[1]) > 0) {
            theta2 = conv * std::acos(1);
        }
    } else {
        theta2 = conv * std::asin(aimpoint_unit_vector[1] / thdenom);
    }

    float vn1 = -1.f * std::cos(theta2 / conv) * std::cos(phi / conv);
    float vn2 = -1.f * std::sin(theta2 / conv) * std::cos(phi / conv);
    float vn3 = std::sin(phi / conv);

    if (aimpoint_unit_vector[0] > 0 && aimpoint_unit_vector[1] > 0 && aimpoint_unit_vector[2] > 0) {
        vn1 = -1.f * std::abs(vn1);
        vn2 = -1.f * std::abs(vn2);
    } else if (aimpoint_unit_vector[0] > 0 && aimpoint_unit_vector[1] > 0 && aimpoint_unit_vector[2] < 0) {
        vn1 = 1.f * std::abs(vn1);
        vn2 = 1.f * std::abs(vn2);
    } else if (aimpoint_unit_vector[0] > 0 && aimpoint_unit_vector[1] < 0 && aimpoint_unit_vector[2] > 0) {
        vn1 = -1.f * std::abs(vn1);
        vn2 = 1.f * std::abs(vn2);
    } else if (aimpoint_unit_vector[0] > 0 && aimpoint_unit_vector[1] < 0 && aimpoint_unit_vector[2] < 0) {
        vn1 = 1.f * std::abs(vn1);
        vn2 = -1.f * std::abs(vn2);
    } else if (aimpoint_unit_vector[0] < 0 && aimpoint_unit_vector[1] > 0 && aimpoint_unit_vector[2] > 0) {
        vn1 = 1.f * std::abs(vn1);
        vn2 = -1.f * std::abs(vn2);
    } else if (aimpoint_unit_vector[0] < 0 && aimpoint_unit_vector[1] > 0 && aimpoint_unit_vector[2] < 0) {
        vn1 = -1.f * std::abs(vn1);
        vn2 = 1.f * std::abs(vn2);
    } else if (aimpoint_unit_vector[0] < 0 && aimpoint_unit_vector[1] < 0 && aimpoint_unit_vector[2] > 0) {
        vn1 = 1.f * std::abs(vn1);
        vn2 = 1.f * std::abs(vn2);
    } else if (aimpoint_unit_vector[0] < 0 && aimpoint_unit_vector[1] < 0 && aimpoint_unit_vector[2] < 0) {
        vn1 = -1.f * std::abs(vn1);
        vn2 = -1.f * std::abs(vn2);
    }

    vn3 = 1.f * std::abs(vn3);
    north_vector[0] = vn1;
    north_vector[1] = vn2;
    north_vector[2] = vn3;

    // calculate the right vector
    right_vector[0] = (aimpoint_unit_vector[1] * vn3) - (aimpoint_unit_vector[2] * vn2);
    right_vector[1] = (aimpoint_unit_vector[2] * vn1) - (aimpoint_unit_vector[0] * vn3);
    right_vector[2] = (aimpoint_unit_vector[0] * vn2) - (aimpoint_unit_vector[1] * vn1);

    plac = (-1.f * dist[0] * position.x) + (-1.f * dist[1] * position.y) + (-1.f * dist[2] * position.z);

    std::vector<float> sun = {23455.f, 0.f, 0.f};
    std::vector<float> p2pout = PointToPlane(sun, plac, 1);
    angle_sun = 270 - (conv * std::atan2(0 - p2pout[0], 0 - p2pout[1]));

    if (angle_sun > 360) {
        angle_sun = angle_sun - 360;
    } else if (angle_sun < 0) {
        angle_sun = angle_sun + 360;
    }
}

void Camera::GetSky()
{
    for (int i = 0; i < nsky; i++) {
        std::vector<float> skygse = Helper::RaDecGSEConversion(std::vector<float> {sky1[i], sky2[i]}, gei_to_gse);
        skygse = Helper::MatrixScalarMultiply(skygse, 1e19);
        std::vector<float> p2pout = PointToPlane(skygse, plac, 1.f);

        if (p2pout[0] == 0.f || p2pout[1] == 0.f) {
            std::cout << "ERROR: RaDecGSEConversion returning invalid values.\n";
        }

        xsky[i] = p2pout[0];
        ysky[i] = p2pout[1];
    }
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
    auto t0 = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> outer = t1 - t0;
    std::chrono::duration<double> xy2radec = t1 - t0;
    std::chrono::duration<double> inner = t1 - t0;
    std::chrono::duration<double> the_rest = t1 - t0;

    auto t_all = std::chrono::high_resolution_clock::now();

    float to_rad = M_PI / 180.f;
    float to_deg = 180.f / M_PI;

    // Let's calculate some angles. Firstly, a rotation about the z axis so that y points towards the aimpoint
    float x_diff = aim.x - position.x;
    float theta = std::atan(x_diff / position.y);
    float rotation_z = ((180 - (to_deg * theta)) * to_rad) + (20 * to_rad);
    float hypotenuse = std::sqrt(std::pow(x_diff, 2) + std::pow(position.y, 2));
    float phi = std::atan(hypotenuse / position.z);
    float rotation_x =  (90 * to_rad) - phi;

    float h = std::sqrt(std::pow(aim.z - position.x, 2) + std::pow(aim.y - position.y, 2));
    float alpha = std::tan(std::abs(h) / std::abs(aim.z - position.z));

    // And now beta
    float beta = (180.f - (std::tan((aim.x - position.x) / position.y) * to_deg)) * to_rad;

    float xh = std::sqrt(std::pow(position.y, 2) + std::pow(position.z, 2));
    float angle = std::atan2(x_diff, xh);

    std::vector<std::vector<float>> rz = {
        {std::cos(rotation_z), - std::sin(rotation_z), 0},
        {std::sin(rotation_z), std::cos(rotation_z), 0},
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

    std::vector<std::vector<float>> rotation = Helper::MatrixMultiply(ry, rx);
    // std::vector<std::vector<float>> rotation2 = Helper::MatrixMultiply(rotation, rx);
    std::vector<std::vector<float>> rotation_inverse = Helper::GetInverse(rotation);

    // Alternative method
    // We first find the angle we have rotated by along the x axis to get to an intermediary frame.
    // Then we calculate the second angle which we have rotated around the z axis.
    // We're using SOHCAHTOA rules here, namely angle = tan(opposite / adjacent)
    // The first angle we will call alpha (a)
    // The second angle will be called beta (b)
    /// Let's calculate alpha
    // float h = std::sqrt(std::pow(aim.z - position.x, 2) + std::pow(aim.y - position.y, 2));
    // float alpha = std::tan((aim.y - position.y) / (aim.z - position.z)) * to_deg;

    // // And now beta
    // float beta = 180.f - std::tan((aim.x - position.x) / position.y) * to_deg;

    // These values remain the same throughout the entire render. These angles only change if the spacecraft position is moved or is the aimpoint has changed.
    // The above calculations also assume that the aimpoint is somewhere along y = 0 and z = 0.
    // We now begin out pixel-specific calculations.

    // serial method first, no parallelisation yet
    for (int i = 0; i < image_dimension; i++) {
        for (int j = 0; j < image_dimension; j++) {

            // Alternative method
            float phi = 180.f + lat[i]; // it is possible that this value is the same for all pixels on the same row as each other...not sure yet

            // The standard formula is x = r cos (phi - theta)
            // where:
            // - r is distance of the line
            // - phi is angle of rotation from line to positive axis
            // - theta is angle between the new and old axis
            // The "axis" is usually x, but in this specific 3d scenario, we are using z. It will be x in the second rotation later on.
            struct {
                float x;
                float y;
                float z;
            } pixel;

            float r = std::sqrt(std::pow(aim.y - position.y, 2) + std::pow(aim.z - position.z, 2));
            pixel.z = r * std::cos((phi - alpha) * to_rad);
            float s = r * std::sin((phi - alpha) * to_rad); // this would be the y coord but it's actualy the magnitude of our next line

            // We will overwrite y, but z is accurate because we will now rotate along the z new z axis, meaning it won't change.
            // Let's calculate our second angle of rtoation here.

            // need to recalculate h

            float theta = 90 - lon[j];
            // float s = std::sqrt(std::pow(aim.x - position.x, 2) + std::pow(aim.y - position.y, 2));
            pixel.x = s * std::cos((theta - beta) * to_rad);
            pixel.y = s * std::sin((theta - beta) * to_rad);


            // Alternative method END

            auto t_outer = std::chrono::high_resolution_clock::now();
            
            auto t_xy2radec = std::chrono::high_resolution_clock::now();
            // Get Right Ascension and Declination for this pixel
            // std::vector<float> pixel_radec = Helper::XYToRaDec(lon[i], lat[j], xsky, ysky, sky1, sky2, skysep);
            auto t_xy2radec_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> d_xy2radec = t_xy2radec_end - t_xy2radec;
            xy2radec = xy2radec + d_xy2radec;

            // Calculate pixel position
            float device_x = (i + 0.5f) / (image_dimension);
            float device_y = (j + 0.5f) / (image_dimension);
            float screen_x = ((2 * device_x) - 1) * std::tan(to_rad * (fov / 2));
            float screen_y = (1 - (2 * device_y)) * std::tan(to_rad * (fov / 2));
            float screen_vec[3] = {screen_x, screen_y, -1};

            float new_distance = Helper::VectorDistance(screen_vec);
            float screen_unit_vec[3] = {
                screen_x / new_distance,
                screen_y / new_distance,
                -1 / new_distance
            }; // That's our camera-based ray unit vector. We need to apply the rotation matrix

            // float vector[3] = {pixel.x, pixel.y, pixel.z};
            
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

            auto t_outer_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> d_outer = t_outer_end - t_outer;
            outer = outer + d_outer;

            for (int pxk = 0; pxk < pxn_dist; pxk++) {
                if (pxk_ok == 1) {

                    auto t_inner = std::chrono::high_resolution_clock::now();
                    
                    // vector to sample point
                    std::vector<float> sample_vector = Helper::MatrixScalarMultiply(world_unit_vector, ray_dist[pxk]);

                    // x coordinate of sample point
                    float x_coord = sample_vector[0] + position.x; // use nearest neighbour
                    float y_coord = sample_vector[1] + position.y;
                    float z_coord = sample_vector[2] + position.z;

                    if ((i == 0 && j == 0) ||
                    (i == 0 && j == 143) ||
                    (i == 143 && j == 0) ||
                    (i == 143 && j == 143)) {
                        if (pxk == 199) {
                            std::cout << pixel.x << "," << pixel.y << "," << pixel.z << "," << std::flush;
                            std::cout << "|||||";
                        }
                    }

                    // if ((x_coord < 8.f && x_coord > 7.5f) &&
                    // (y_coord < .25f && y_coord > -0.25f) &&
                    // (z_coord < .25f && z_coord > -0.25f)
                    // ) {
                    //     image[i][j] = 7;
                    //     continue;
                    // }

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
                    sample = sample * ray_width * 1000;

                    image[i][j] = image[i][j] + sample;
                    sample_vector.clear();

                    pxk_yes = 1;

                    auto t_inner_end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> d_inner = t_inner_end - t_inner;
                    inner = inner + d_inner;
                }
            }
            // pixel_radec.clear();
            // pixel_ray_vector.clear();
            // unit_vector.clear();
            
        }
        std::cout << std::to_string(i) << ", " << std::flush;
    }

    auto t_all_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d_all = t_all_end - t_all;

    std::cout << "\nImage complete.\n";
    std::cout << "Outer loop:" << std::to_string(outer.count()) << "\n";
    std::cout << "XY2RaDec:" << std::to_string(xy2radec.count()) << "\n";
    std::cout << "Inner loop:" << std::to_string(inner.count()) << "\n";
    std::cout << "Rendering:" << std::to_string(d_all.count()) << "\n";
    

}

int Camera::ToFITS()
{
    long naxis = 2;
    long naxes[2] = {144, 144};

    // using namespace CCfits;

    // const CCfits::String& filename = "testfits";
    // CCfits::RWmode mode = CCfits::RWmode::Write;

    // CCfits::FITS pFits = CCfits::FITS(filename, 100, naxis, naxes);

    std::ofstream outfile("../python/testfile.dat");

    for (int i = 0; i < 144; i++) {
        for (int j = 0; j < 144; j++) {
            outfile << std::to_string(image[i][j]);
            outfile << ",";
        }
        outfile << "\n";
    }

    return 1;
}

Camera::~Camera()
{
}
