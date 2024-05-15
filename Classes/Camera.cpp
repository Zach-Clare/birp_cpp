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
    // TODO: GetRayStepDistance(33)
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
    float apon = conv * std::cos(t3); // not sure what apon means (angle in plane of north?)

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
    int soy = 1;
    int sox = 1;
    if (apor > 90) {
        sox = -1;
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
        t3 = 1;
    }
    float observer_object_angle = conv * std::acos(t3);

    float xo = 0.f;
    float yo = 0.f;
    if (observer_object_angle > 0.000001f) {
        xo = sox * observer_object_angle * std::sin(apon / conv);
        yo = soy * observer_object_angle * std::cos(apon / conv);
    }

    return {xo, yo, aro, observer_object_angle};
}

void Camera::Orient()
{
    float conv = 180 / M_PI;

    dist[0] = aim.x - position.x;
    dist[1] = aim.y - position.y;
    dist[2] = aim.z - position.z;

    float distance = sqrt((pow(dist[0], 2)) + (pow(dist[1], 2)) + (pow(dist[2], 2)));

    // aimpoint vector as a unit vector
    aimpoint_unit_vector[0] = dist[0] / distance;
    aimpoint_unit_vector[1] = dist[1] / distance;
    aimpoint_unit_vector[1] = dist[2] / distance;
    
    // calculate the north vector
    float phi = conv * std::acos(aimpoint_unit_vector[2]);
    float thdenom = std::sin(phi / conv);

    // calculate based around unit vectors x value
    float theta1;
    if (abs(thdenom) < abs(aimpoint_unit_vector[0])) {
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
    if (abs(thdenom) < abs(aimpoint_unit_vector[1])) {
        if ((thdenom * aimpoint_unit_vector[1]) < 0) {
            theta2 = conv * std::acos(-1);
        } else if ((thdenom * aimpoint_unit_vector[1]) > 0) {
            theta2 = conv * std::acos(1);
        }
    } else {
        theta2 = conv * std::acos(aimpoint_unit_vector[1] / thdenom);
    }

    float vn1 = -1.f * std::cos(theta2 / conv) * std::cos(phi / conv);
    float vn2 = -1.f * std::sin(theta2 / conv) * std::cos(phi / conv);
    if (aimpoint_unit_vector[0] > 0 && aimpoint_unit_vector[1] > 0 && aimpoint_unit_vector[2] > 0) {
        vn1 = -1.f * abs(vn1);
        vn2 = -1.f * abs(vn2);
    } else if (aimpoint_unit_vector[0] > 0 && aimpoint_unit_vector[1] > 0 && aimpoint_unit_vector[2] < 0) {
        vn1 = 1.f * abs(vn1);
        vn2 = 1.f * abs(vn2);
    } else if (aimpoint_unit_vector[0] > 0 && aimpoint_unit_vector[1] < 0 && aimpoint_unit_vector[2] > 0) {
        vn1 = -1.f * abs(vn1);
        vn2 = 1.f * abs(vn2);
    } else if (aimpoint_unit_vector[0] > 0 && aimpoint_unit_vector[1] < 0 && aimpoint_unit_vector[2] < 0) {
        vn1 = 1.f * abs(vn1);
        vn2 = -1.f * abs(vn2);
    } else if (aimpoint_unit_vector[0] < 0 && aimpoint_unit_vector[1] > 0 && aimpoint_unit_vector[2] > 0) {
        vn1 = 1.f * abs(vn1);
        vn2 = -1.f * abs(vn2);
    } else if (aimpoint_unit_vector[0] < 0 && aimpoint_unit_vector[1] > 0 && aimpoint_unit_vector[2] < 0) {
        vn1 = -1.f * abs(vn1);
        vn2 = 1.f * abs(vn2);
    } else if (aimpoint_unit_vector[0] < 0 && aimpoint_unit_vector[1] < 0 && aimpoint_unit_vector[2] > 0) {
        vn1 = 1.f * abs(vn1);
        vn2 = 1.f * abs(vn2);
    } else if (aimpoint_unit_vector[0] < 0 && aimpoint_unit_vector[1] < 0 && aimpoint_unit_vector[2] < 0) {
        vn1 = -1.f * abs(vn1);
        vn2 = -1.f * abs(vn2);
    }

    float vn3 = 1.f * abs(vn3);
    north_vector[0] = vn1;
    north_vector[1] = vn2;
    north_vector[2] = vn3;

    // calculate the right vector
    right_vector[0] = (aimpoint_unit_vector[1] * vn3) - (aimpoint_unit_vector[2] * vn2);
    right_vector[1] = (aimpoint_unit_vector[2] * vn1) - (aimpoint_unit_vector[0] * vn1);
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
        std::vector<float> p2pout = PointToPlane(skygse, plac, 1);

        if (p2pout[0] == 0.f || p2pout[1] == 0.f) {
            std::cout << "ERROR: RaDecGSEConversion returning invalid values.";
        }

        xsky[i] = p2pout[0];
        ysky[i] = p2pout[1];
    }
}

Camera::~Camera()
{
}
