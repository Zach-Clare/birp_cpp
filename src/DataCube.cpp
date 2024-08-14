#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include "Helper.h"

#include "DataCube.h"

void DataCube::Load(std::string filename, int data_begin, bool debug = false) {
	std::cout << "Loading datacube...\n";

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

	std::cout << "Datacube loaded.\n";
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