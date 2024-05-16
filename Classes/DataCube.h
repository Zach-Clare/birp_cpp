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
	struct{
		int x;
		int y;
		int z;
	} size;
	
	float spacing[3];
	float origin_coords[3];

	std::vector<std::vector<std::vector<float>>> slices; // needs to be of length k (or coords_z?)

private:
	void InitSpacing();
	void InitOriginCoords();
	void InitSize();
	float GetAxisSpacing(std::vector<float> coords);
};
