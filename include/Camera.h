/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


#pragma once
#include "DataCube.h"

//!  Camera can generate sightlines then sample from the given Space object 

/*!
  Camera aims to simulate the SXI camera located onboard the SMILE spacecraft.
  It will generate rays and smaple the given Space object, whether that's a
  DataCube-like object or a CMEM-like function.
*/
class Camera
{
private:
    int image_dimension; //!< in pixels, how tall and wide the resultant render should be
    int image_dimension_x;
    int image_dimension_y;
    float fov_x; //!< Used with pixel_size_deg to calculate fov of th render
    float fov_y;
    struct {
        float x;
        float y;
        float z;
    } position; //!< Position of the camera in GSE space
    struct {
        float x;
        float y;
        float z;
    } aim; //!< Aimpoint of the centre pixel. Only supports non-zero values for x. Non-zero y and z values will be ignored.
    std::vector<float> ray_dist; //!< Vector of ray distances to multiply with the ray unit vectors to get sample point
    std::vector<float> ray_widths; //!< If enabled, can modulate the width of the rays to account for the spready of the rays as they diverge from the camera
    float pixel_size_deg; //!< Along with FOV, decides the field of view of the camera
    float pixel_size_rad; //!< Calculated from pizel_size_width
    float ray_width;
    int ray_samples; //!< Number of samples to take along each ray

    void GenerateRayDistWidth(int);
    void Integrate();
    float Orient();
    std::vector<float> PointToPlane(std::vector<float>, float, float, std::vector<float> distance_vec, std::vector<float> unit_vec, std::vector<float> north, std::vector<float> right);

public:
    Space* dataCube; // dataCube is a misleading name since it could be a CMEM object
    // float image[144][144] {0}; // image is square, it's size will be image_dimension^2 
    std::vector<std::vector<float>> image;

    // Note cube is passed by reference for efficiency
    Camera(Space& cube, float pixel_size_deg, float plot_fov_h, float plot_fov_w);
    ~Camera();

    void SetPosition(float&, float&, float&); //!< Set the position of the camera in GSE space
    void SetAim(float, float, float); //!< Set the aimpoint of the camera in GSE space. Only non-zero x values are supported
    void Render(); //!< Begin the rendering process
    int ToDat(std::string filename); //!< Export to a simple CSV file with a .dat file extension
    // int ToFITS(std::string filename); //!< Eport the image to a FITS file
};