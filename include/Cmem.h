#pragma once
#include "Space.h"

class CMEM : public Space {
public:
    void Init();
    float GetSample(float x, float y, float z);

private:
    float LinScaled(
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
    );
    static std::vector<float> ShueCoords(float x, float y, float z);

};