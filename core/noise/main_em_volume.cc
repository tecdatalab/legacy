// Load several mrc files and return the volume for each file.
// All voxels

#include "density_lattice.h"
#include "mrc.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv) 
{
    if(argc < 3) {
		cerr << "For each file, the program returns the volume corresponding to the voxels that have "
			 << "a density higher or equal than density_threshold" << endl << endl;
        cerr << "Usage: " << argv[0] << " density_threshold mrc_file1 mrc_file2 ..." << endl;
        return 0;
    }
	int current_file_index = 2;
	float density_threshold = atof(argv[1]);
	while (current_file_index < argc) {
	    mrc_reader* reader = new mrc_reader(argv[current_file_index]);
    	density_lattice* lattice = reader->get_densities();
		double volume = lattice->calculate_volume(density_threshold);
    float min, max;
    lattice->get_min_max_densities(&min, &max);
		cout << argv[current_file_index] << " "
         << std::fixed << volume;
    cout.unsetf ( std::ios::floatfield );
    cout << " nx: " << lattice->size_x()
         << " ny: " << lattice->size_y()
         << " nz: " << lattice->size_z() << "\n\n"
         << "Min/Max densities: " << min << " " << max << "\n\n"
         << "voxel sizes: "
         << lattice->voxel_length_x() << " "
         << lattice->voxel_length_y() << " "
         << lattice->voxel_length_z() << "\n\n"
         << "origin: "
         << reader->get_x_origin() << " "
         << reader->get_y_origin() << " "
         << reader->get_z_origin() << endl;

    delete lattice;
   	delete reader;
		current_file_index++;
	}
	return 0;
}
