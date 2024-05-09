#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include "../helper.cpp"

class DataCube {
public:

	std::vector<float> coords_x;
	std::vector<float> coords_y;
	std::vector<float> coords_z;
	
	float spacing[3];

	void Load(std::string filename, int data_begin) {
		std::cout << "Loading data...\n";

		std::fstream ifs(filename); // create filestream from filename
		std::string line; // create space for line data to go
		int line_num = 0; // start at line zero
		std::string meta;
		std::vector<float> line_buffer;
		std::vector<float> coord_buffer;
		int buffer_i = 0;

		while(std::getline(ifs, line)) {
			++line_num; // incrememnt line number

			if (line_num == 1) {
				meta = line;
			}

			// this area describes the coord array
			if ((11 <= line_num) && (line_num <= (data_begin - 1))) {
				line_buffer = explode_float(line, ' ');
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

			// coord arrays built. Now build slices

			if (line_num == data_begin) {
				this->GetSpacingAll();
			}
			int thing = 1+1;
			
		}
		std::cout << std::to_string(this->spacing[0]);
		std::cout << "\ndeeeeeef\n";
	}

private:

	void GetSpacingAll() {
		this->spacing[0] = this->GetSpacing(coords_x);
		this->spacing[1] = this->GetSpacing(coords_y);
		this->spacing[2] = this->GetSpacing(coords_z);
	}

	float GetSpacing(std::vector<float> coords) {
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
		if (!std::equal(distances.begin() + 1, distances.end(), distances.begin())) {
			//all equal
			throw std::invalid_argument("non-uniform grid given");
		}

		return distances[0];
	}
};