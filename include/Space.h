#pragma once

class Space {
public:
    // pure virtual functions for abstract class
    virtual void Init() = 0;
    virtual float GetSample(float x, float y, float z) = 0;

};