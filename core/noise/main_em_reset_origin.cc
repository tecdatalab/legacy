#include <cstdlib>
#include <iostream>
#include <vector>
#include "../proteinstructure/atom.h"
#include "density_lattice.h"
#include "main_em_reset_origin_options.h"
#include "mrc.h"
#include "../proteinstructure/pdb.h"

using std::cerr;
using std::endl;
using std::string;

int main(int argc, char** argv)
{
	main_em_reset_origin_options options(argc,argv);
	if(!options.parse_successful())
	{
    cout << "Adjusts the X,Y,Z origin fields from EMAN generated MRC files.\n"
         << "The adjustment is meant to allow the exact alignment of the \n"
         << "density map and the PDB coordinates that originated it.\n\n";
    cout << "Usage: em_reset_origin --simulated-em <mrc file> --base-pdb "
         << "<pdb file> --contour-level <surface isovalue> --output <mrc file>"
         << endl;
    cout << "       em_reset_origin --simulated-em <mrc file> "
         << "-x <number> -y <number> -z <number> --output <mrc file>"
         << endl;
		exit(EXIT_FAILURE);
	}

	string em_file = options.get_em_file();
  string pdb_file = options.get_pdb_file();
	string output_file = options.get_output_file();
  float contour_level = options.get_contour_level();

	// Load the information
  mrc_reader* simulated_em_reader = new mrc_reader(em_file);
  density_lattice* simulated_em_lattice = simulated_em_reader->get_densities();

  if(options.is_direct_origin_provided()) {
    // then it's assumed no PDB is used and we just write a new file with
    // the explicit values provided
    cout << "Writing to " << output_file << "...\n";

    simulated_em_reader->set_origin(options.get_x(), options.get_y(),
                                    options.get_z());
    simulated_em_reader->write(output_file);
 
    delete simulated_em_lattice;
    delete simulated_em_reader;

	  exit(EXIT_SUCCESS);
  }
  // Else there's not explicit new origin, infer from the PDB coordinates.
	pdb input_pdb;
	read_protein(pdb_file, input_pdb);
  // Get the boundaries in voxel space, for the EM map, and in PDB coordinate
  // space, for the PDB file.
  density_voxel min_em_bound, max_em_bound;
  float min_pdb_x, min_pdb_y, min_pdb_z, max_pdb_x, max_pdb_y, max_pdb_z;

  simulated_em_lattice->get_voxel_boundaries_for_threshold(
      contour_level, &min_em_bound, &max_em_bound);
  input_pdb.get_boundaries(&min_pdb_x, &min_pdb_y, &min_pdb_z,
                           &max_pdb_x, &max_pdb_y, &max_pdb_z);

  // EMAN first translates the lowest x,y,z to 0,0,0 and then adds empty
  // voxel padding around the volume.
  // Set the origin in the output EM to compensate for these 2 adjustments.
  unsigned long voxel_padding_x =
      /* total voxels - (voxels in/on the surface boundary) */
      simulated_em_lattice->size_x() - (max_em_bound.x - min_em_bound.x);
  unsigned long voxel_padding_y =
      simulated_em_lattice->size_y() - (max_em_bound.y - min_em_bound.y);
  unsigned long voxel_padding_z =
      simulated_em_lattice->size_z() - (max_em_bound.z - min_em_bound.z);

  // Padding divided by two because padding is added on both sides of each axis
  float padding_adjustment_x = static_cast<float>(voxel_padding_x)
                               * simulated_em_lattice->voxel_length_x() / 2;
  float padding_adjustment_y = static_cast<float>(voxel_padding_y)
                               * simulated_em_lattice->voxel_length_y() / 2;
  float padding_adjustment_z = static_cast<float>(voxel_padding_z)
                               * simulated_em_lattice->voxel_length_z() / 2;

  // The adjustment up to this point has been in real space, but in the file
  // the origin is set in terms of the number of voxels, so we need to divide
  // by the length of each axis.
  float adjustment_x = (min_pdb_x - padding_adjustment_x) /
                       simulated_em_lattice->voxel_length_x();
  float adjustment_y = (min_pdb_y - padding_adjustment_y) /
                       simulated_em_lattice->voxel_length_x();
  float adjustment_z = (min_pdb_z - padding_adjustment_z) /
                       simulated_em_lattice->voxel_length_x();

  cout << "Adjustments (X,Y,Z)" << endl
       << "Padding: (" << padding_adjustment_x << ", "
       << padding_adjustment_y << ", " << padding_adjustment_z << ")" << endl
       << "PDB origin translation: (" << min_pdb_x << ", " << min_pdb_y
       << ", " << min_pdb_z << ")" << endl
       << "Total: (" << adjustment_x << ", " << adjustment_y << ", "
       << adjustment_z << ")" << endl;

  cout << "Writing to " << output_file << "...\n";

  simulated_em_reader->set_origin(adjustment_x, adjustment_y, adjustment_z);
  simulated_em_reader->write(output_file);
 
  delete simulated_em_lattice;
  delete simulated_em_reader;

	exit(EXIT_SUCCESS);
}
