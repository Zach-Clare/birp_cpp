#include <vector>
#include <cmath>

#include "DataCube.h"
#include "../Helper.h"

#include "Camera.h"

Camera::Camera(DataCube cube, float pixel_size_deg, int plot_fov)
{
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
        sky1.push_back(c_lon + skysep);

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
}

std::vector<float> Camera::PointToPlane(std::vector<float> obj, float plac, float prad)
{
    float conv = 180 / M_PI;

    float snumer = (dist[0] * obj[0]) + (dist[1] * obj[1]) + (dist[2] * obj[2]) + plac;
    float sdenom = (dist[0] * aimpointUnitVector[0]) + (dist[1] * aimpointUnitVector[1]) + (dist[2] * aimpointUnitVector[2]);
    float scalar = snumer / sdenom;

    float q[3] = {
        obj[0] - (scalar * aimpointUnitVector[0]),
        obj[1] - (scalar * aimpointUnitVector[1]),
        obj[2] - (scalar * aimpointUnitVector[2])
    };

    // vector of object from the camera
    float object_vector[3] = {
        q[0] - position.x,
        q[1] - position.y,
        q[2] - position.z
    };

    float dp = std::sqrt((object_vector[0] * object_vector[0]) + (object_vector[1] * object_vector[1]) + (object_vector[2] * object_vector[2]));
    float object_unitvector[3] = {
        object_vector[0] / dp,
        object_vector[1] / dp,
        object_vector[2] / dp
    };
    dp = std::sqrt((object_unitvector[0] * object_unitvector[0]) + (object_unitvector[1] * object_unitvector[1]) + (object_unitvector[2] * object_unitvector[2]));

    // maps to line 530 Camera.py
}

void Camera::Orient()
{
    float conv = 180 / M_PI;

    dist[0] = aim.x - position.x;
    dist[1] = aim.y - position.y;
    dist[2] = aim.z - position.z;

    float distance = sqrt((pow(dist[0], 2)) + (pow(dist[1], 2)) + (pow(dist[2], 2)));

    // aimpoint vector as a unit vector
    aimpointUnitVector[0] = dist[0] / distance;
    aimpointUnitVector[1] = dist[1] / distance;
    aimpointUnitVector[1] = dist[2] / distance;
    
    // calculate the north vector
    float phi = conv * std::acos(aimpointUnitVector[2]);
    float thdenom = std::sin(phi / conv);

    // calculate based around unit vectors x value
    float theta1;
    if (abs(thdenom) < abs(aimpointUnitVector[0])) {
        if ((thdenom * aimpointUnitVector[0]) < 0) {
            theta1 = conv * std::acos(-1);
        } else if ((thdenom * aimpointUnitVector[0]) > 0) {
            theta1 = conv * std::acos(1);
        }
    } else {
        theta1 = conv * std::acos(aimpointUnitVector[0] / thdenom);
    }

    // calculate based around unit vectors y value
    float theta2;
    if (abs(thdenom) < abs(aimpointUnitVector[1])) {
        if ((thdenom * aimpointUnitVector[1]) < 0) {
            theta2 = conv * std::acos(-1);
        } else if ((thdenom * aimpointUnitVector[1]) > 0) {
            theta2 = conv * std::acos(1);
        }
    } else {
        theta2 = conv * std::acos(aimpointUnitVector[1] / thdenom);
    }

    float vn1 = -1.f * std::cos(theta2 / conv) * std::cos(phi / conv);
    float vn2 = -1.f * std::sin(theta2 / conv) * std::cos(phi / conv);
    if (aimpointUnitVector[0] > 0 && aimpointUnitVector[1] > 0 && aimpointUnitVector[2] > 0) {
        vn1 = -1.f * abs(vn1);
        vn2 = -1.f * abs(vn2);
    } else if (aimpointUnitVector[0] > 0 && aimpointUnitVector[1] > 0 && aimpointUnitVector[2] < 0) {
        vn1 = 1.f * abs(vn1);
        vn2 = 1.f * abs(vn2);
    } else if (aimpointUnitVector[0] > 0 && aimpointUnitVector[1] < 0 && aimpointUnitVector[2] > 0) {
        vn1 = -1.f * abs(vn1);
        vn2 = 1.f * abs(vn2);
    } else if (aimpointUnitVector[0] > 0 && aimpointUnitVector[1] < 0 && aimpointUnitVector[2] < 0) {
        vn1 = 1.f * abs(vn1);
        vn2 = -1.f * abs(vn2);
    } else if (aimpointUnitVector[0] < 0 && aimpointUnitVector[1] > 0 && aimpointUnitVector[2] > 0) {
        vn1 = 1.f * abs(vn1);
        vn2 = -1.f * abs(vn2);
    } else if (aimpointUnitVector[0] < 0 && aimpointUnitVector[1] > 0 && aimpointUnitVector[2] < 0) {
        vn1 = -1.f * abs(vn1);
        vn2 = 1.f * abs(vn2);
    } else if (aimpointUnitVector[0] < 0 && aimpointUnitVector[1] < 0 && aimpointUnitVector[2] > 0) {
        vn1 = 1.f * abs(vn1);
        vn2 = 1.f * abs(vn2);
    } else if (aimpointUnitVector[0] < 0 && aimpointUnitVector[1] < 0 && aimpointUnitVector[2] < 0) {
        vn1 = -1.f * abs(vn1);
        vn2 = -1.f * abs(vn2);
    }

    float vn3 = 1.f * abs(vn3);
    float vn[3] = {vn1, vn2, vn3};

    // calculate the right vector
    float vr[3] = {
        (aimpointUnitVector[1] * vn3) - (aimpointUnitVector[2] * vn2),
        (aimpointUnitVector[2] * vn1) - (aimpointUnitVector[0] * vn1),
        (aimpointUnitVector[0] * vn2) - (aimpointUnitVector[1] * vn1)
    };

    float plac = (-1.f * dist[0] * position.x) + (-1.f * dist[1] * position.y) + (-1.f * dist[2] * position.z);

    float sun[3] = {23455.f, 0.f, 0.f};


}

Camera::~Camera()
{
}
