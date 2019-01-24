/*
 * main_test_noise_removal.cc
 *
 * Simple main to test noise removal functions.
 */

#include "density_lattice.h"
#include "mrc.h"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <input map> <contour isovalue>\n\n"
         << "Noise removal result output: <input map>.mod" << endl;
    return 0;
  }
	//Read in EMD file
	cout << "Reading in file " << argv[1] << "\n";
	mrc_reader* reader = new mrc_reader(argv[1]);     //input file

	//Set densities over suggested contour level from EMD
	cout << "Setting densities " << atof(argv[2]) << "\n";
   	density_lattice* lattice = reader->get_densities(atof(argv[2])); // suggested contour level for 1010

    //Write file back after noise has been removed
    cout << "Writing back after noise removal\n";
    reader->write(argv[1] + string(".mod"), lattice);   //output file

    //Clean up
    delete lattice;
   	delete reader;

	return 0;
}





