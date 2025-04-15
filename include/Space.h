/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


#pragma once

//! An interface defining what the Camera class requires to "look" through the magnetopause.

/*! This is a pure virtual function, meaning it's an interface for the DataCube and
CMEM class to implement. It's very simple, and the Camera class needs to know nothing
more about what it's looking at, other than the functions defiend here. */

class Space {
public:
    // pure virtual functions for abstract class
    virtual void Init() = 0;
    virtual float GetSample(float x, float y, float z) = 0;

};