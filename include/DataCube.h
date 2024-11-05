#pragma once
#include <string>
#include <iostream>
#include <vector>

#include "Space.h"

class DataCube : public Space {
public:
	void Init();
	void Init(std::string& filename, int data_start, bool debug);
	void Load(std::string filename, int data_start, bool debug); //!< Parse an exported MHD file into memory, ready to be sampled. 
	float GetSample(float x, float y, float z); //!< Returns the value from the cube at x, y, and z. Set interpolation with SetTrilinear().
	void SetTrilinear(bool); //!< Tells the sampling function to use Trilinear. If set to false, use nearest neighbour. 

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
	bool trilinear;
};
