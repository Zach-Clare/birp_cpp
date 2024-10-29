#include "Space.h"
#include "Cmem.h"
#include "Helper.h"

#include <iostream>
#include <cmath>
#include <vector>

/* This will initialise function parameters necessary for the CMEM calculations.
*/
void CMEM::Init() {

    // These numbers are magic numbers from Sam's code. I don't know what these numbers do, but I know their value. Or at least Sam does.

    // v is some sort of vector that we need to calculate the dynamic pressure
    v[0] = 400; // x
    v[1] = 0;   // y
    v[2] = 0;   // z

    // z is v but for magnetic pressure
    b[0] = 0;   // x
    b[1] = 0;   // y
    b[2] = 5;   // z

    CalcDynamicPressure();
    CalcMagneticPressure();

    dipole = 0.f;

    // Initialise parameters for...models? Probbaly Shue or Jorgensen?
    p0 = 1;
    bs = r0_lin + 3;
    A1 = 0.000015f;
    A2 = 0.000013f;
    B = 2.f;
    alpha = 2.5f;
    beta = -1.6f;
    p1 = 1;
    p2 = 3;
    p3 = 4;

    // New init methods
    bs = -0.12 * density + 13.24;
    A1 = 0.0000027 * density - 0.0000063;
    A2 = 0.0000009 * density - 0.0000010;
    p0 = 0.0022 * density + 0.7753;
    
    CalcInitialAlpha();

    DefineLinearCoeffs();

    return;
}

float CMEM::GetSample(float x, float y, float z) {

    // convert {x, y, z} to {r, theta, and phi}
    std::vector<float> shue = CMEM::ShueCoords(x, y, z);

    float radius_mp = LinScaled(shue[1], shue[2], dn, ds, theta_n, theta_s, r0_lin, p0, p1, p2, p3);
    float radius_bs = ShueModel(shue[1], shue[2], bs, ay_bs, az_bs);

    if (shue[0] < radius_mp) {
        return 0;
    } else if (shue[0] >= radius_mp && (shue[0] < radius_bs)) {
        return A1 * (std::exp(- B * std::pow(shue[1] / 2, 4))) * std::pow(shue[0] / 10, (- alpha - (beta * std::pow(std::sin(shue[1]), 2))));
    } else {
        return A2 * (std::pow(shue[0] / 10, -3));
    }
}

float CMEM::LinScaled(
        float& theta,
        float& phi,
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

        float phi_n = std::acos(
            (std::cos(theta) * std::cos(theta_n)) +
            (std::sin(theta) * std::sin(theta_n) * std::cos(phi - M_PI / 2))
        );
        float phi_s = std::acos(
            (std::cos(theta) * std::cos(theta_s)) +
            (std::sin(theta) * std::sin(theta_s) * std::cos(phi - (3 * M_PI / 2)))
        );

        // what does f mean??
        float f = std::pow(std::cos(theta / 2) + a[5] * std::sin(2 * theta) * (1 - std::exp( - theta)), (p1 * (beta_c[0] + beta_c[1] * std::cos(phi) + beta_c[3] * std::pow(std::sin(phi), 2))));

        // what is q?
        float q = p2 * charlie * std::exp(dn * (std::pow(phi_n, a[21]))) + p2 * charlie * std::exp(ds * (std::pow(phi_s, a[21])));

        float r = p0 * r0_lin * f + q;

        return r;
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

/* Defines Linear Coefficients stored in CMEM's private member variables
*/
void CMEM::DefineLinearCoeffs() {

    // alpha coeefficients
    a = {12.544f, -0.194f, 0.305f, 0.0573f, 2.178f, 0.0571f, -0.999f, 16.473f, 0.00152f, 0.382f, 0.0431f, -0.00763f, -0.210f, 0.0405f, -4.430f, -0.636f, -2.600f, 0.832f, -5.328f, 1.103f, -0.907f, 1.450f};

    // from these and func args, get beta coefficients
    beta_c = {
        a[6] + a[7] * ((std::exp(a[8] * b[2]) - 1) / (std::exp(a[9] * b[2]) + 1)),
        a[10],
        a[11] + a[12] * dipole, a[13]
    };

    // next, get charlie n and s coeffs, which are equal
    charlie = a[14] * std::pow((pressure_dynamic + pressure_magnetic), a[15]);

    // get delta coeffs - north and south?
    // not identical, look carefully
    dn = a[16] + a[17] * dipole + a[18] * std::pow(dipole, 2); // note +
    ds = a[16] - a[17] * dipole + a[18] * std::pow(dipole, 2); // note -

    // get theta coeffs - north and south?
    theta_n = a[19] + a[20] * dipole; // note +
    theta_s = a[19] - a[20] * dipole; // note -

    // get unscaled subsolar magnetopause radius
    r0_lin = a[0] * std::pow(pressure_dynamic + pressure_magnetic, a[1]) * (1 + a[2] * ((std::exp(a[3] * b[2]) - 1) / (std::exp(a[4] * b[2]) + 1)));
}

/* Implementation of the Shue model defined in Jorgensen et al. (2019). Don't ask me about the maths.
*/
float CMEM::ShueModel(float theta, float phi, float r0, float ay, float az) {
    float ry = r0 * std::pow((2 / (1 + std::cos(theta))), ay);
    float rz = r0 * std::pow((2 / (1 + std::cos(theta))), az);

    float r = (ry * rz) / std::pow(std::pow(rz * std::cos(phi), 2) + std::pow(ry * std::sin(phi), 2), 0.5);
    return r;
}

/* Calculated value is constant throughout simulation (I think)
*/
void CMEM::CalcDynamicPressure() {
    // assuming average ion mass = mass of proton
    float proton_mass = 0.00000000000000000000000000167f;

    // calculate v ? in m/s
    float vms = Helper::VectorDistance(v) * 1000;

    // Convert number of particles from cm^3 to m^3
    float n = density * 1000000; // That's 1,000,000 - a million

    // Calculate dynamic pressure first in Pascals, then nPa
    pressure_dynamic = 0.5f * proton_mass * n * std::pow(vms, 2) * 1000000000; // That's 1,000,000,000 - a billion
}

/* Calculated value is constant throughout simulation (I think)
*/
void CMEM::CalcMagneticPressure() {
    // Calculate magnitude of B in T ?
    B = Helper::VectorDistance(b) * 0.000000001;

    // ?
    float mu0 = 4 * M_PI * 0.0000001;

    // Calculate magnetic pressue in Pascals, then nPa
    pressure_magnetic = std::pow(B, 2) / (2 * mu0) * 1000000000; // That's 1,000,000,000 - a billion
    
}

void CMEM::CalcInitialAlpha() {
    float val = (0.58f - 0.010 * b[2]) * (1 + 0.010f * pressure_dynamic) + 0.2f;

    ay_bs = val;
    az_bs = val;
}