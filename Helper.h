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
        static std::vector<std::vector<float>> Gei2GseTransforms();
        static std::vector<float> RaDecToGSI(float ra, float dec);
        static std::vector<float> RaDecGSEConversion(std::vector<float>, std::vector<std::vector<float>>);
        
};