#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>

#include "Helper.h"

#include "DataCube.h"

void DataCube::Init() {
	return;
}

void DataCube::Init(std::string &filename, int data_begin, bool debug = false) {
	this->Load(filename, data_begin, debug);
}

void DataCube::Load(std::string filename, int data_begin, bool debug) {

	std::fstream ifs(filename); // create filestream from filename
	std::string line; // create space for line data to go
	int line_num = 0; // start at line zero
	std::string meta;
	std::vector<float> line_buffer;
	std::vector<float> coord_buffer;
	int buffer_i = 0;
	std::vector<float> data_buffer;

	int i = 0;
	int j = 0;
	int k = 0;

	slices.push_back(std::vector<std::vector<float>> {});

	while(std::getline(ifs, line)) {
		++line_num; // incrememnt line number

		if (line_num == 1) {
			meta = line;
		}

		// this area describes the coord array
		if ((11 <= line_num) && (line_num <= (data_begin - 1))) {
			line_buffer = Helper::explode_float(line, ' ');
			coord_buffer.insert(coord_buffer.end(), line_buffer.begin(), line_buffer.end());

			if (line_buffer.size() < 6) {
				if (buffer_i == 0) {
					coords_x = coord_buffer;
				} else if (buffer_i == 1)
				{
					coords_y = coord_buffer;
				} else if (buffer_i == 2)
				{
					coords_z = coord_buffer;
				}

				++buffer_i;
				coord_buffer.clear();
			} else {
				line_buffer.clear();
			}
		}

		// we have all coord values. Initialise meta-coord values
		if (line_num == data_begin) {
			this->InitSpacing();
			this->InitOriginCoords();
			this->InitSize();
		}

		// All coord-related arrays built. Now build slices with actual data
		if (line_num >= (data_begin - 1)) {
			// trim line. 
			Helper::trim(line);
			if (line == "") {
				k++;
				j = 0;
				slices.push_back(std::vector<std::vector<float>> {});
				continue;
			}

			// process volume
			if (line_num >= data_begin) {
				std::vector<float> inter = Helper::explode_float(line, '  '); // split line into individual elements
				data_buffer.insert(data_buffer.end(), inter.begin(), inter.end());

				if (inter.size() < 6) { // end of this thread, meaning data buffer is full. Add the thread
					slices[k].push_back(data_buffer);
				
					data_buffer.clear();
					j++;
				}
			}
		}
		
	}

	if (debug == true) {
		slices[80][79][50] = 1;

		for (int i = 0; i < 101; i++) {
			slices[80][79][i] = 1;
		}
	}

	if (slices.size() == 0) {
		throw std::invalid_argument("Datacube is invalid, inaccessible, or malformed.");
	}
}

void DataCube::InitSpacing() {
	this->spacing[0] = this->GetAxisSpacing(coords_x);
	this->spacing[1] = this->GetAxisSpacing(coords_y);
	this->spacing[2] = this->GetAxisSpacing(coords_z);
}

float DataCube::GetAxisSpacing(std::vector<float> coords) {
	// first calculate the distances between each element
	auto size = coords.size(); // get the size of coords
	std::vector<float> distances; // create the holding array
	for (int i = 0; i < size; i++) { // loop through each of coords items
		if (i == 0) { // first element, skip
			continue;
		}

		// use abs() because sometimes we deal with negative coordinates
		distances.push_back(std::abs(coords[i - 1] - coords[i]));
	}

	// now check that all the distances are the same
	// shamelessly taken from Raxvan via StackOverflow:
	// https://stackoverflow.com/questions/20287095/checking-if-all-elements-of-a-vector-are-equal-in-c
	if (!std::equal(distances.begin() + 1, distances.end(), distances.begin())) {
		//all equal
		throw std::invalid_argument("non-uniform grid given");
	}

	return distances[0];
}

void DataCube::InitOriginCoords() {
	// Simply fill out the array with the origin values
	// (does this not need to be the coords of where the zero point are?)
	origin_coords[0] = this->coords_x[0];
	origin_coords[1] = this->coords_y[0];
	origin_coords[2] = this->coords_z[0];
}

void DataCube::InitSize() {
	size.x = this->coords_x.size();
	size.y = this->coords_y.size();
	size.z = this->coords_z.size();
}

void DataCube::SetTrilinear(bool option) {
	trilinear = option;
}

float DataCube::GetSample(float x, float y, float z) {
	// If this gives a segmentation fault, the datacube likely has a problem loading. TODO: Error catch
	float x_ingress = std::abs(coords_x[0] - x);
	float x_index = (x_ingress / spacing[0]);

	float y_ingress = std::abs(coords_y[0] - y);
	float y_index = (y_ingress / spacing[1]);

	float z_ingress = std::abs(coords_z[0] - z);
	float z_index = (z_ingress / spacing[2]);

	if ((0 > x_index || x_index + 1 >= size.x) ||
		(0 > y_index || y_index + 1 >= size.y) ||
		(0 > z_index || z_index + 1 >= size.z)) {
		return -1.f;
	}
	
	float sample;
	if (trilinear) {
		// we need the 8 values of our 1x1x1 cube to interpolate through
		// so we'll floor then ceiling our float index.
		int z_min = std::floor(z_index);
		int z_max = std::ceil(z_index);
		int y_min = std::floor(y_index);
		int y_max = std::ceil(y_index);
		int x_min = std::floor(x_index);
		int x_max = std::ceil(x_index);

		// Don't be scared, defining the points of the cube is just
		// counting in binary!
		
		// Define the various points
		std::vector<int> p000 = {z_min, y_min, x_min};
		std::vector<int> p001 = {z_min, y_min, x_max};
		std::vector<int> p010 = {z_min, y_max, x_min};
		std::vector<int> p011 = {z_min, y_max, x_max};
		std::vector<int> p100 = {z_max, y_min, x_min};
		std::vector<int> p101 = {z_max, y_min, x_max};
		std::vector<int> p110 = {z_max, y_max, x_min};
		std::vector<int> p111 = {z_max, y_max, x_max};

		// Retrieve the values at each of these points
		float s000 = slices.at(p000[0]).at(p000[1]).at(p000[2]);
		float s001 = slices.at(p001[0]).at(p001[1]).at(p001[2]);
		float s010 = slices.at(p010[0]).at(p010[1]).at(p010[2]);
		float s011 = slices.at(p011[0]).at(p011[1]).at(p011[2]);
		float s100 = slices.at(p100[0]).at(p100[1]).at(p100[2]);
		float s101 = slices.at(p101[0]).at(p101[1]).at(p101[2]);
		float s110 = slices.at(p110[0]).at(p110[1]).at(p110[2]);
		float s111 = slices.at(p111[0]).at(p111[1]).at(p111[2]);

		// get the fractional part of the float (e.g. 4.352 -> 0.352)
		float x_frac = x_index - (int) x_index;
		float y_frac = y_index - (int) y_index;
		float z_frac = z_index - (int) z_index;

		// Interpolate between first dimension x (four times)
		float low_x1 = std::lerp(s000, s001, x_frac);
		float low_x2 = std::lerp(s010, s011, x_frac);
		float high_x1 = std::lerp(s100, s101, x_frac);
		float high_x2 = std::lerp(s110, s111, x_frac);

		// Interpolate second dimension y (twice)
		float low_y = std::lerp(low_x1, low_x2, y_frac);
		float high_y = std::lerp(high_x1, high_x2, y_frac);

		// Interpolate final third dimension z (once)
		sample = std::lerp(low_y, high_y, z_frac);

	} else {
		sample = slices.at((int) z_index).at((int) y_index).at((int) x_index);
	}
	
	return sample;
}