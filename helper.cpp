#include <string>
#include <vector>
#include <sstream>
#include <utility>

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