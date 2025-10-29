/*
-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------
*/


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

void DataCube::Init(std::string &filename, bool debug = false) {
	this->Load(filename, debug);
}

void DataCube::Load(std::string filename, bool debug) {

	std::fstream ifs(filename); // create filestream from filename
	std::string line; // create space for line data to go
	int line_num = 0; // start at line zero
	std::string meta;
	std::vector<float> line_buffer; // This holds values from each line
	std::vector<float> coord_buffer; // This holds all values belonging to one axis (and gets reset)
	std::vector<float> data_buffer;
	std::vector<int> coord_array;

	int i = 0;
	int j = 0;
	int k = 0;

	slices.push_back(std::vector<std::vector<float>> {});

	bool parsed_coords = false;
	bool check_parse = false;
	bool parsed_space = false;

	while(std::getline(ifs, line)) {
		++line_num; // incrememnt line number

		if (line_num == 1) {
			meta = line;
		}

		if (line_num == 10){
			coord_array = Helper::explode_int(line, ' ');
			check_parse = true; // this prevents a segfault when we check against the intended coord size

		}

		// this area describes the coord array
		if (check_parse // check parse prevents segfault
			&& coords_x.size() == coord_array[0]
			&& coords_y.size() == coord_array[1]
			&& coords_z.size() == coord_array[2]) {
			parsed_coords = true; // use this as a marker. Previously, this used a line number, but it was very inelastic
		}
		
		// now we have the intended size of the datacube, start parsing the coordinate values
		if ((11 <= line_num) && !parsed_coords) {
			line_buffer = Helper::explode_float(line, ' '); // get this line and...
			coord_buffer.insert(coord_buffer.end(), line_buffer.begin(), line_buffer.end()); // add it to the overall buffer

			if (coords_x.size() == 0 && coord_buffer.size() == coord_array[0]) { // if no x coords yet and coord buffer is correct size, use it.
				coords_x = coord_buffer;
				coord_buffer.clear();
			} else if (coords_y.size() == 0 && coord_buffer.size() == coord_array[1]) { // or maybe it's time to fill Y?
				coords_y = coord_buffer;
				coord_buffer.clear();
			} else if (coords_z.size() == 0 && coord_buffer.size() == coord_array[2]) { // or z?
				coords_z = coord_buffer;
				coord_buffer.clear();
			} else {
				// No action required, just keep building coord_buffer
			}
		}

		// we have all coord values. Initialise meta-coord values
		// if we do this without checking if we've parsed the coordinates, we'll try and parse the space before we've built any of the necessary arrays.
		if (parsed_coords && !parsed_space) {
			this->InitSpacing();
			this->InitOriginCoords();
			this->InitSize();
			parsed_space = true;
		}

		// All coord-related arrays built. Now build slices with actual data
		if (parsed_coords) { // make sure we're at the correct part of the input file - all coordinate arrays should be built by now
			Helper::trim(line); // trim line with custom function // TODO: REmove????
			
			if (slices[k].size() == coord_array[1]) { // if our slices are the correct size, move  onto filling the next slice
				k++; // This increases in z direction I think
				j = 0; // reset Y counter
				slices.push_back(std::vector<std::vector<float>> {}); // create new vector for next set of values
				if (line == "") {
					continue;
				}
			}

			// process volume
			if (parsed_coords) {
				std::vector<float> inter = Helper::explode_float(line, '  '); // split line into individual elements
				data_buffer.insert(data_buffer.end(), inter.begin(), inter.end());

				if (data_buffer.size() == coord_array[0]) { // end of this thread, meaning data buffer is full. Add the thread
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
	// if (!std::equal(distances.begin() + 1, distances.end(), distances.begin())) {
	if (!Helper::EqualDistance(distances)) {
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