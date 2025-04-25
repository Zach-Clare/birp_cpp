/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


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
#include <cfloat>

#include "Helper.h"

std::vector<float> Helper::explode_float(std::string const& s, char delim)
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

std::vector<int> Helper::explode_int(std::string const& s, char delim)
{
    std::vector<int> result;
    std::istringstream iss(s);

    for (std::string token; std::getline(iss, token, delim); )
    {
        if (!token.empty()) {
            result.push_back(std::move(std::stoi(token)));
        }
    }

    return result;
}

std::vector<std::string> explode_string(std::string const & s, char delim)
{
    std::vector<std::string> result;
    std::stringstream iss(s);

    for (std::string token; std::getline(iss, token, delim);)
    {
        if (!token.empty()) {
            result.push_back(std::move(token));
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
    // int columns = b.size(); // b is 1D
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

std::vector<float> Helper::ApplyRotation(std::vector<std::vector<float>> rot, float* coords)
{
    std::vector<float> result;
    // we want to use coords vertically, as if there are three nested single-element arrays
    // First, loop over rot
    for (auto row : rot) {
        // take each element and multiply it by the equivalent in coords
        // and add it to a buffer
        result.push_back((row[0] * coords[0]) + (row[1] * coords[1]) + (row[2] * coords[2]));
    }

    return result;
}

float Helper::GetDeterminant(const std::vector<std::vector<float>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } 
    int dimension = vect.size();

    if(dimension == 0) {
        return 1;
    }

    if(dimension == 1) {
        return vect[0][0];
    }

    //Formula for 2x2-matrix
    if(dimension == 2) {
        return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
    }

    float result = 0;
    int sign = 1;
    for(int i = 0; i < dimension; i++) {

        //Submatrix
        std::vector<std::vector<float>> subVect(dimension - 1, std::vector<float> (dimension - 1));
        for(int m = 1; m < dimension; m++) {
            int z = 0;
            for(int n = 0; n < dimension; n++) {
                if(n != i) {
                    subVect[m-1][z] = vect[m][n];
                    z++;
                }
            }
        }

        //recursive call
        result = result + sign * vect[0][i] * Helper::GetDeterminant(subVect);
        sign = -sign;
    }

    return result;
}

std::vector<std::vector<float>> Helper::GetTranspose(const std::vector<std::vector<float>> matrix1) {

    //Transpose-matrix: height = width(matrix), width = height(matrix)
    std::vector<std::vector<float>> solution(matrix1[0].size(), std::vector<float> (matrix1.size()));

    //Filling solution-matrix
    for(size_t i = 0; i < matrix1.size(); i++) {
        for(size_t j = 0; j < matrix1[0].size(); j++) {
            solution[j][i] = matrix1[i][j];
        }
    }
    return solution;
}

std::vector<std::vector<float>> Helper::GetCofactor(const std::vector<std::vector<float>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } 

    std::vector<std::vector<float>> solution(vect.size(), std::vector<float> (vect.size()));
    std::vector<std::vector<float>> subVect(vect.size() - 1, std::vector<float> (vect.size() - 1));

    for(std::size_t i = 0; i < vect.size(); i++) {
        for(std::size_t j = 0; j < vect[0].size(); j++) {

            int p = 0;
            for(size_t x = 0; x < vect.size(); x++) {
                if(x == i) {
                    continue;
                }
                int q = 0;

                for(size_t y = 0; y < vect.size(); y++) {
                    if(y == j) {
                        continue;
                    }

                    subVect[p][q] = vect[x][y];
                    q++;
                }
                p++;
            }
            solution[i][j] = pow(-1, i + j) * Helper::GetDeterminant(subVect);
        }
    }
    return solution;
}

std::vector<std::vector<float>> Helper::GetInverse(const std::vector<std::vector<float>> vect) {
    if(Helper::GetDeterminant(vect) == 0) {
        throw std::runtime_error("Determinant is 0");
    } 

    float d = 1.f / Helper::GetDeterminant(vect);
    std::vector<std::vector<float>> solution(vect.size(), std::vector<float> (vect.size()));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            solution[i][j] = vect[i][j]; 
        }
    }

    solution = Helper::GetTranspose(Helper::GetCofactor(solution));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            solution[i][j] *= d;
        }
    }

    return solution;
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
    // using namespace std::chrono_literals;

    float conv = M_PI / 180; // conversion from degrees to radians

    float date = DateToJulianDate(2024, 1, 1);
    date = date - 2400000.5f; // modify to be relative to 17th Nov 1858
    float t0 = (date - 51544.5f) / 36525; // convert into julian centuries

    // float theta = 100.461f + 36000.770f * t0;
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

std::vector<float> Helper::TwoMin(std::vector<float> arr, std::vector<int> mask) 
{
    float smallest = 0;
    float second = 1;
    if (arr[second] < arr[smallest]) {
        float temp = smallest;
        smallest = second;
        second = temp;
    }
    for (int i = 2; i < arr.size() - 1; i++) {
        if (std::find(mask.begin(), mask.end(), i) == mask.end()) {
            continue; // don't search this element
        }
        if (arr[i] < arr[smallest]) {
            second = smallest;
            smallest = i;
        } else if (arr[i]  < arr[second]) {
            second = i;
        }
    }

    return std::vector<float> {smallest, second};
}

std::vector<float> Helper::XYToRaDec(float ct_x, float ct_y, std::vector<float> xsky, std::vector<float> ysky, std::vector<float> sky1, std::vector<float> sky2, float skysep)
{
    float conv = 180 / M_PI;
    float out[2] = {-999.999f, -999.999f};

    // float ra = std::atan2(ct_x, ct_y);

    std::vector<float> dtest;
    for (int i = 0; i < xsky.size(); i++) { 
        dtest.push_back(std::sqrt(((ct_x - xsky[i]) * (ct_x - xsky[i])) + ((ct_y - ysky[i]) * (ct_y - ysky[i]))));
    } // calculate closest value in xsky/ysky to ct_x/ct_y
    int dmax = std::distance(dtest.begin(), std::min_element(dtest.begin(), dtest.end()));

    std::vector<int> sel;
    for (int i = 0; i < sky1.size(); i++) {
        if (sky1[i] == sky1[dmax]) {
            sel.push_back(i);
        }
    }

    // std::vector<float> sel_values;
    // for (int i = 0; i < sel.size(); i++) {
    //     sel_values.push_back(dtest[sel[i]]);
    // }

    std::vector<float> min_two = TwoMin(dtest, sel);
    float diff_dec = sky2[min_two[1]] - sky2[min_two[0]];
    float diff_x = xsky[min_two[1]] - xsky[min_two[0]];
    float diff_y = ysky[min_two[1]] - ysky[min_two[0]];

    sel.clear();
    for (int i = 0; i < sky2.size(); i++) {
        if (sky2[i] == sky2[dmax]) {
            sel.push_back(i);
        }
    }
    min_two = TwoMin(dtest, sel);
    float diff_ra = sky1[min_two[1]] - sky1[min_two[0]];
    float diff_x2 = xsky[min_two[1]] - xsky[min_two[0]];
    float diff_y2 = ysky[min_two[1]] - ysky[min_two[0]];

    float d_ang = conv * std::atan2(diff_y, diff_x);
    float c_ang = conv * std::atan2((ct_y - ysky[dmax]), (ct_x - xsky[dmax]));
    float ct_dec = sky2[dmax] + ((diff_dec / skysep) * (dtest[dmax] * std::cos((d_ang - c_ang) / conv)));

    if (ct_dec > 90.f) {
        ct_dec = 90.f;
    } else if (ct_dec < -90.f) {
        ct_dec = -90.f;
    }

    if (diff_ra == (skysep - 360.f)) {
        diff_ra = skysep;
    } else if (diff_ra == (360.f - skysep)) {
        diff_ra = -1.f * skysep;
    }
    float d2_ang = conv * std::atan2(diff_y2, diff_x2);

    float ct_ra;
    if (std::cos(ct_dec / conv) != 0) {
        ct_ra = sky1[dmax] + ((diff_ra / (1.f * skysep)) * (dtest[dmax] * std::cos((d2_ang - c_ang) / conv)) / std::cos(ct_dec / conv));
    } else {
        ct_ra = sky1[dmax];
    }

    if (ct_ra < 0) {
        ct_ra = ct_ra + 360.f;
    } else if (ct_ra > 360) {
        ct_ra = ct_ra - 360.f;
    }

    return std::vector<float> {ct_ra, ct_dec};
}

std::vector<float> Helper::Interp1( std::vector< float > &x, std::vector< float > &y, std::vector< float > x_new )
{
    std::vector< float > y_new;
    y_new.reserve( x_new.size() );

    std::vector< float > dx, dy, slope, intercept;
    dx.reserve( x.size() );
    dy.reserve( x.size() );
    slope.reserve( x.size() );
    intercept.reserve( x.size() );
    for( int i = 0; i < x.size(); ++i ){
        if( i < x.size()-1 )
        {
            dx.push_back( x[i+1] - x[i] );
            dy.push_back( y[i+1] - y[i] );
            slope.push_back( dy[i] / dx[i] );
            intercept.push_back( y[i] - x[i] * slope[i] );
        }
        else
        {
            dx.push_back( dx[i-1] );
            dy.push_back( dy[i-1] );
            slope.push_back( slope[i-1] );
            intercept.push_back( intercept[i-1] );
        }
    }

    for ( int i = 0; i < x_new.size(); ++i ) 
    {
        int idx = findNearestNeighbourIndex( x_new[i], x );
        y_new.push_back( slope[idx] * x_new[i] + intercept[idx] );

    }

    return y_new;

}

int Helper::findNearestNeighbourIndex( float value, std::vector< float > &x )
{
    float dist = FLT_MAX;
    int idx = -1;
    for ( int i = 0; i < x.size(); ++i ) {
        float newDist = value - x[i];
        if ( newDist > 0 && newDist < dist ) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

std::vector<float> Helper::XYZToSpherical(float* xyz)
{
    float x = xyz[0];
    float y = xyz[1];
    float z = xyz[2];

    float r = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));

    float theta;
    // if (z > 0.f) {
        theta = std::atan2(std::sqrt(std::pow(x, 2) + std::pow(y, 2)), z);
    // } else if (z == 0.f && std::sqrt(std::pow(x, 2) + std::pow(y, 2)) != 0.f) {
    //     theta = M_PI / 2;
    // }

    float phi;
    phi = std::atan2(y, x);

    return {phi, theta, r};
}

std::vector<float> Helper::XYZToSphericalAlt(float* xyz)
{
    float x = xyz[0];
    float y = xyz[1];
    float z = xyz[2];


    float r = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));

    float theta;
    // if (z > 0) {
    //     theta = std::atan(std::sqrt(std::pow(x, 2) + std::pow(y, 2)) / z);
    // } else if (z < 0) {
    //     theta = M_PI + std::atan(std::sqrt(std::pow(x, 2) + std::pow(y, 2)) / z);
    // } else if (z == 0.f && std::sqrt(std::pow(x, 2) + std::pow(y, 2)) != 0) {
    //     theta = M_PI / 2;
    // }
    // theta = std::atan2(std::sqrt(std::pow(x, 2) + std::pow(y, 2)), z);
    theta = std::atan2(y, x);

    float phi;
    // if (x > 0) {
    //     phi = std::atan(y / x);
    // } else if (x < 0 && y >= 0) {
    //     phi = std::atan(y / x) + M_PI;
    // } else if (x < 0 && y < 0) {
    //     phi = std::atan(y / x) - M_PI;
    // } else if (x == 0.f && y > 0) {
    //     phi = M_PI / 2;
    // } else if (x == 0.f && y < 0) {
    //     phi = -(M_PI / 2);
    // }
    // phi = std::acos(z / r);
    phi = std::atan2(z, std::sqrt(std::pow(x, 2) + std::pow(y, 2)));

    return {phi, theta, r};
}

// Calculate distance of 3D vector
float Helper::VectorDistance(float* xyz)
{
    return std::sqrt(std::pow(xyz[0], 2) + std::pow(xyz[1], 2) + std::pow(xyz[2], 2));
}

// Calculate distance of 2D vector
float Helper::VectorDistance2D(float* xy)
{
    return std::sqrt(std::pow(xy[0], 2) + std::pow(xy[1], 2));
}

// steps should be a list of distanecs from one element of the list/vector to the next
bool Helper::EqualDistance(std::vector<float> steps)
{
    // we need to take this list of distances and account for floating point errors.
    // to do this, we can loop through each element and see if it's within a certain distance of it's neighbour.
    for (std::vector<float>::iterator i = steps.begin(); i != steps.end(); ++i) {
        // as soon as one of them isn't we return false
        if (i == steps.end() - 1) {
            // final element
            break;
        }

        // continue;

        if (i - std::next(i) > 0.00001f) {
            return false;
        }
    }
    
    // if that doesn't happen, we return true
    return true;
}