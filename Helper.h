#pragma once
#include <string>
#include <vector>
#include <algorithm>

// Trim functions work in place
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

        static std::vector<float> explode_float(std::string const & s, char delim);
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
        static std::vector<float> XYZToSpherical(float*, int);
        static std::vector<float> XYZToSphericalAlt(float*, int);
        static float VectorDistance(float* xyz);
        
};