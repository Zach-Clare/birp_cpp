#pragma once
#include <string>
#include <iostream>
#include <vector>

class DataCube {
public:
	void Load(std::string, int);

	std::vector<float> coords_x;
	std::vector<float> coords_y;
	std::vector<float> coords_z;
	
	float spacing[3];
	float origin_coords[3];

private:
	void InitSpacing();
	void InitOriginCoords();
	float GetAxisSpacing(std::vector<float> coords);
};
