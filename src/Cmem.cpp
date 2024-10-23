#include "Space.h"
#include "Cmem.h"
#include "Helper.h"

#include <iostream>
#include <cmath>

void CMEM::Init() {
    std::cout << "please" << std::endl; // please :(
    return;
}

float CMEM::GetSample(float x, float y, float z) {

    // convert {x, y, z} to {r, theta, and phi}
    std::vector<float> shue = CMEM::ShueCoords(x, y, z);


    return 0.1f; // filler
}

float CMEM::LinScaled(
        float& theta,
        float& phi,
        float& a,
        float& beta_c,
        float& c,
        float& dn,
        float& ds,
        float& theta_n,
        float& theta_s,
        float& r0_lin,
        float& p0,
        float& p1,
        float& p2,
        float& p3
    ) {
        // This function calculates r using the Lin model.
}

std::vector<float> CMEM::ShueCoords(float x, float y, float z) {
    float point[3] = {x, y, z};
    float phi_components[2] = {y, z};

    float r = Helper::VectorDistance(point);

    float theta;
    float phi;

    if (x != 0) {
        theta =  std::acos(x / r);
    } else {
        theta = 0;
    }

    if (y != 0) {
        phi = std::acos(y / Helper::VectorDistance2D(phi_components));
    } else {
        phi = 0;
    }

    return std::vector<float> {r, theta, phi};
}
