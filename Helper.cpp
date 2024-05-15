#include <string>
#include <vector>
#include <sstream>
#include <utility>
#include <algorithm>
#include <cctype>
#include <locale>
#include <math.h>
#include <chrono>
#include <ctime>
#include <iterator>

#include "Helper.h"

std::vector<float> Helper::explode_float(std::string const & s, char delim)
{
    std::vector<float> result;
    std::istringstream iss(s);

    for (std::string token; std::getline(iss, token, delim); )
    {
        if (!token.empty()) {
            result.push_back(std::move(std::stof(token)));
        }
    }

    return result;
}

std::vector<std::vector<float>> Helper::MatrixMultiply(std::vector<std::vector<float>> a, std::vector<std::vector<float>> b)
{
    // WARNING
    // assuming regular (rows = columns) matraces
    // and assuming symmetrical matraces
    int size_a = a.size();
    int size_b = b.size();

    std::vector<std::vector<float>> result;

    for (int i = 0; i < size_a; i++) {
        std::vector<float> buffer;
        for (int j = 0; j < size_b; j++) {
            // result.push_back({0, 0, 0, 0});
            float value = 0.0f;

            for (int k = 0; k < size_a; k++) {
                // result[i][j] += a[i][k] * b[k][j];
                value = value + (a[i][k] * b[k][j]);
            }

            buffer.push_back(value);
        }
        result.push_back(buffer);
        buffer.clear();
    }

    return result;
}

std::vector<float> Helper::MatrixMultiply(std::vector<std::vector<float>> a, std::vector<float> b) {
    int rows = a.size(); // a is 2D
    int a_columns = a[0].size();
    int columns = b.size(); // b is 1D
    float result[4] = {0}; // fixed size
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < a_columns; col++) {
            result[row] += a[row][col] * b[col];
        }
    }

    std::vector<float> v(std::begin(result), std::end(result));

    return v;
}

std::vector<std::vector<float>> Helper::MatrixScalarMultiply(std::vector<std::vector<float>> matrix, float scalar)
{
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[0].size(); j++) {
            matrix[i][j] = matrix[i][j] * scalar;
        }
    }

    return matrix;
}

std::vector<float> Helper::MatrixScalarMultiply(std::vector<float> matrix, float scalar)
{
    for (int i = 0; i < matrix.size(); i++) {
        matrix[i] = matrix[i] * scalar;
    }

    return matrix;
}

float DateToJulianDate(int year, int month, float day) 
{
    int yearp;
    int monthp;
    if (month == 1 || month == 2){
        yearp = year - 1;
        monthp = month + 12;
    } else {
        yearp = year;
        monthp = month;
    }

    // check where we are in relation to October 15 1582, which
    // is when we switched from Julian to Gregorian
    int a;
    int b;
    int c;
    if ((year < 1582) ||
        (year == 1582 && month < 10) ||
        (year == 1582 && month == 10 && day < 15)) {
            b = 0;
    } else {
        a = std::trunc(yearp / 100.f);
        b = 2 - a + std::trunc(a / 4.f);
    }

    if (yearp < 0) {
        c = std::trunc((356.25 * yearp) - 0.75);
    } else {
        c = std::trunc(365.25 * yearp);
    }

    int d = std::trunc(30.6001 * (monthp + 1));

    float jd = b + c + d + day + 1720994.5f;
    return jd;
}

std::vector<std::vector<float>> Helper::Gei2GseTransforms()
{
    using namespace std::chrono_literals;

    float conv = M_PI / 180; // conversion from degrees to radians

    float date = DateToJulianDate(2024, 1, 1);
    date = date - 2400000.5f; // modify to be relative to 17th Nov 1858
    float t0 = (date - 51544.5f) / 36525; // convert into julian centuries

    // auto now = std::chrono::system_clock::now(); // time right now
    // auto julian_date = now + 58574100h; // shift to julian calendar
    // auto julian_date_mod = julian_date - (2400000.5f * 24h); // modified relative to 17 Nov 1858
    // auto point = std::chrono::time_point_cast<std::chrono::hours>(julian_date_mod); // Fancy footwork. cast to time_point
    // auto epoch = point.time_since_epoch(); // get duration from time_point object
    // auto duration = std::chrono::duration_cast<std::chrono::hours>(epoch); // cast to duration
    // float t0 = (epoch - (52544.5f * 24h)).count() / 36525; // so we can use .count() to get time in julian centuries

    float theta = 100.461f + 36000.770f * t0;
    float epsilon = 23.439f - 0.013f * t0;
    float m = 357.528 + 35999.050f * t0;
    float lambda = 280.460 + 36000.772 * t0;
    float lambda_sun = lambda + (1.915 - 0.0048 * t0) * std::sin(m * conv) + 0.020f * std::sin((2 * m) * conv);

    std::vector<std::vector<float>> gei2gseA = {
        {std::cos(lambda_sun * conv), std::sin(lambda_sun * conv), 0, 0},
        {-std::sin(lambda_sun * conv), std::cos(lambda_sun * conv), 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };

    std::vector<std::vector<float>> gei2gseB {
        {1, 0, 0, 0},
        {0, std::cos(epsilon * conv), std::sin(epsilon * conv), 0},
        {0, -std::sin(epsilon * conv), std::cos(epsilon * conv), 0},
        {0, 0, 0, 1}
    };

    return MatrixMultiply(gei2gseA, gei2gseB);
}

std::vector<float> Helper::RaDecToGSI(float ra, float dec)
{
    float conv = M_PI / 180; // conversion float

    float alpha = ra * conv;
    float delta = dec * conv;

    float z_gei = std::sin(delta);
    float y_int = 1 - std::pow(std::sin(delta), 2);
    if (y_int < 0) {
        y_int = 0.f;
    }

    float y_gei = std::sqrt(y_int) * std::tan(alpha) / std::sqrt(1 + std::pow(std::tan(alpha),2));
    float x_int = 1 - std::pow(y_gei, 2) - std::pow(std::sin(delta), 2);
    if (x_int < 0) {
        x_int = 0.f;
    }

    float x_gei = std::sqrt(x_int);

    if (90 < ra && ra < 180) {
        x_gei = std::abs(x_gei) * (-1);
        y_gei = std::abs(y_gei);
    } else if (180 < ra && ra < 270) {
        x_gei = std::abs(x_gei) * (-1);
        y_gei = std::abs(y_gei) * (-1);
    }

    return std::vector<float> {x_gei, y_gei, z_gei};
}

std::vector<float> Helper::RaDecGSEConversion(std::vector<float> ra_dec, std::vector<std::vector<float>> gei_to_gse)
{
    std::vector<float> gei_co = RaDecToGSI(ra_dec[0], ra_dec[1]);
    gei_co.push_back(1.f);
    std::vector<float> gse_pos = MatrixMultiply(gei_to_gse, gei_co);
    gse_pos.pop_back();

    return gse_pos;
}