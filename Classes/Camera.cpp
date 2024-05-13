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

    std::cout << "Sky generation complete\n";

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

Camera::~Camera()
{
}
