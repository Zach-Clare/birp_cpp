#pragma once
#include "Space.h"

#include <vector>

class CMEM : public Space {
public:
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
        float& charlie,
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
    void CalcInitialAlpha();

};