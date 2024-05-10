#pragma once
#include <string>
#include <vector>
#include <sstream>
#include <utility>
#include <algorithm>
#include <cctype>
#include <locale>

// Trim functions work in place
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

inline void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}

std::vector<float> explode_float(std::string const & s, char delim)
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