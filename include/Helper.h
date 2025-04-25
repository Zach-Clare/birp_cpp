/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


#pragma once
#include <string>
#include <vector>
#include <algorithm>

//! A collection of functions that don't belong to any of the existing classes.

/*! Think of this as a workbook, or a collection of generic math helpers. */

class Helper 
{
    public:
        // Trim functions work in place
        inline static void ltrim(std::string &s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
        }

        inline static void rtrim(std::string &s) {
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), s.end());
        }

        inline static void trim(std::string &s) {
            rtrim(s);
            ltrim(s);
        }

        // ends_with courtesy of https://stackoverflow.com/a/2072890
        inline static bool ends_with(std::string const & value, std::string const & ending) {
            if (ending.size() > value.size()) return false;
            return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
        }

        static std::vector<float> explode_float(std::string const & s, char delim);
        static std::vector<int> explode_int(std::string const & s, char delim);
        static std::vector<std::string> explode_string(std::string const & s, char delim);
        static std::vector<std::vector<float>> MatrixMultiply(std::vector<std::vector<float>> a, std::vector<std::vector<float>> b);
        static std::vector<float> MatrixMultiply(std::vector<std::vector<float>> a, std::vector<float> b);
        static std::vector<std::vector<float>> MatrixScalarMultiply(std::vector<std::vector<float>>, float);
        static std::vector<float> MatrixScalarMultiply(std::vector<float>, float);
        static std::vector<float> ApplyRotation(std::vector<std::vector<float>> rot, float* coords);
        static std::vector<std::vector<float>> GetInverse(const std::vector<std::vector<float>> vect);
        static float GetDeterminant(const std::vector<std::vector<float>> vect);
        static std::vector<std::vector<float>> GetTranspose(const std::vector<std::vector<float>> matrix1);
        static std::vector<std::vector<float>> GetCofactor(const std::vector<std::vector<float>> vect);
        static std::vector<std::vector<float>> Gei2GseTransforms();
        static std::vector<float> RaDecToGSI(float ra, float dec);
        static std::vector<float> RaDecGSEConversion(std::vector<float>, std::vector<std::vector<float>>);
        static std::vector<float> TwoMin(std::vector<float>, std::vector<int>);
        static std::vector<float> XYToRaDec(float, float, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, float);
        static std::vector<float> Interp1( std::vector< float > &x, std::vector< float > &y, std::vector< float > x_new );
        static int findNearestNeighbourIndex( float value, std::vector< float > &x );
        static std::vector<float> XYZToSpherical(float*);
        static std::vector<float> XYZToSphericalAlt(float*);
        static float VectorDistance(float* xyz);
        static float VectorDistance2D(float* xy);
        static bool EqualDistance(std::vector<float> steps);
        
};