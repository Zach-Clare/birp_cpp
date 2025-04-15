/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


#pragma once
#include "Space.h"

#include <vector>

//! CMEM is a sample-able function space

/*! An alternative to a datacube, CMEM generates emission values mathematically.*/

class CMEM : public Space {
public:
    void Init(
        bool hack,
        std::vector<int> v_passed = {}, 
        std::vector<int> b_passed = {}, 
        float dipole_passed = 0.f,
        float p0_passed = 0.786300004f,
        int p1_passed = 1.f,
        int p2_passed = 3.f,
        int p3_passed = 4.f,
        float B_passed = 2.f,
        float alpha_passed = 2.5f,
        float beta_passed = -1.6f,
        float bs_passed = 12.6400003f,
        float A1_passed = 7.2e-06f,
        float A2_passed = 3.5000000000000004e-06f,
        float ay_bs_passed = NULL,
        float az_bs_passed = NULL);
    void Init();
    float GetSample(float x, float y, float z);

private:

    // Input parameters (I think?)
    int temp = 200000; // That's 200,000
    float density = 5;
    float v[3];
    float b[3];
    float pressure_dynamic;
    float pressure_magnetic;
    float dipole;

    // Linear Coefficients
    std::vector<float> a;
    std::vector<float> beta_c;
    float charlie;
    float dn;
    float ds;
    float theta_n;
    float theta_s;
    float r0_lin;

    // Other values ?
    float p0; // scaling factor on the subsolar magnetopause parameter
    float bs; // subsolar bowshock distance parameter
    float A1; // parameter ?
    float A2; // parameter ?
    float B; // parameter ?
    float alpha; // parameter ?
    float beta; // parameter ?
    float p1; // scalaing factor for magnetopause flaring parameter
    float p2; // scaling factor for magnetopause indentation parameter
    float p3; // scaling factor for ? parameter
    float ay_bs; // ay bowshock flaring parameter
    float az_bs; // az bowshock flaring parameter

    float LinScaled(
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
    );
    static std::vector<float> ShueCoords(float x, float y, float z);
    void DefineLinearCoeffs();
    static float ShueModel(float theta, float phi, float r0, float ay, float az);
    void CalcDynamicPressure();
    void CalcMagneticPressure();
    float CalcInitialAlpha();

};