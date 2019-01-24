#include <cstdlib>
#include <iostream>
#include <vector>
#include "align.h"
#include "atom.h"
#include "main_pdb_align_options.h"
#include "pdb.h"
#include "point_transformation.h"

using std::cerr;
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
	main_pdb_align_options options(argc,argv);
	if(!options.parse_successful())
	{
    cerr << "Aligns --input atomic coordinates to the matching atoms in "
         << "--align_to, according to the residue name, chain id and sequence "
         << "identifier, and outputs the modified coordinates to a new file. "
         << "It only uses c-alpha atoms for the alignment." << endl;
    cerr << "It is expected that all atoms in --input have a match in --align-to." << endl;
    cerr << "If --transformed-pdb is provided the atoms transformed and "
         << "output come from this file, not from --input" << endl;
    cerr << "If --ignore-chain-id is provided only residue names and residue id"
         << " values are used for matching purposes."  << endl << endl;
		cerr << "Usage: pdb_align --input <pdb file> --align-to <pdb file> --output <pdb file> [--transformed-pdb <pdbfile>] [--ignore-chain-id]" << endl;
		exit(EXIT_FAILURE);
	}

	string input_file = options.get_input_file();
  string output_file = options.get_output_file();
	string align_to_file = options.get_align_to_file();
  string transformed_file = options.get_transformed_file();
  bool ignore_chain_id = options.get_ignore_chain_id();

	// Load the information
	pdb input_pdb, align_to_pdb;

	read_protein(input_file, input_pdb);
	read_protein(align_to_file, align_to_pdb);

  vector<atom> input_matches;
  vector<atom> align_to_matches;

  input_pdb.get_matching_calphas(&align_to_pdb, ignore_chain_id, &input_matches,
                                 &align_to_matches);

  if(input_matches.size() != align_to_matches.size() ||
     (!input_matches.size())) {
    cerr << "Error: Incompatible alignment. All atoms in --input must match "
         << "an atom in --align-to" << endl;
		exit(EXIT_FAILURE);
  }

  vector<atom> transformed_atoms;

  // Check if the optional transformed-pdb parameter was provided
  // and use the atoms in that file, if that's the case
  if (!transformed_file.empty()) {
    pdb transformed_pdb;
    read_protein(transformed_file, transformed_pdb);
    transformed_atoms = transformed_pdb.atoms;
  } else {
    transformed_atoms = input_pdb.atoms;
  }

  vector<atom> aligned_output;
  point_transformation optimal_transformation;
  double rmsd = align_all_and_transform(
      &input_matches, &align_to_matches, &transformed_atoms,
      &aligned_output, &optimal_transformation);

  cout << "C-alpha atoms aligned: "<< input_matches.size() << "\nRMSD: " << rmsd << endl;
  cout << "Writing to " << output_file << endl;

  write_complex(output_file, aligned_output);

#ifdef DEBUG
/*  for(size_t i =0; i < input_matches.size(); i++) {
    atom matched_atom = align_to_matches[i];
    atom current_calpha = input_matches[i];
    atom aligned_atom = aligned_output[i];
        cout << "XXX " << matched_atom.atype << " " << matched_atom.rnum << " " << matched_atom.residue
             << " " << matched_atom.axyz[0] << " " << matched_atom.axyz[1] << " " << matched_atom.axyz[2] << endl;
        cout << "YYY " << current_calpha.atype << " " << current_calpha.rnum << " " << current_calpha.residue
             << " " << current_calpha.axyz[0] << " " << current_calpha.axyz[1] << " " << current_calpha.axyz[2] << endl;
        cout << "ZZZ " << aligned_atom.atype << " " << aligned_atom.rnum << " " << aligned_atom.residue
             << " " << aligned_atom.axyz[0] << " " << aligned_atom.axyz[1] << " " << aligned_atom.axyz[2] << endl;
  }
*/
  cout << "Input atoms " << input_pdb.atoms.size() << " Aligned Output " << aligned_output.size() << endl;
  double x,y,z,rx,ry,rz;
  for(size_t i =0; i < input_pdb.atoms.size(); i++) {
    atom input_atom = input_pdb.atoms[i];
    atom aligned_atom = aligned_output[i];
    optimal_transformation.transform(input_atom.axyz[0], input_atom.axyz[1], input_atom.axyz[2], &x, &y, &z);
    optimal_transformation.invert_transform(x,y,z,&rx,&ry,&rz);
    cout << aligned_atom.atype << " " << aligned_atom.rnum << " " << aligned_atom.residue
         << " " << aligned_atom.axyz[0] << " " << aligned_atom.axyz[1] << " " << aligned_atom.axyz[2] << " --- "
         << x << " " << y << " " << z << " +++ "
         << " " << input_atom.axyz[0] << " " << input_atom.axyz[1] << " " << input_atom.axyz[2] << " --- "
         << rx << " " << ry << " " << rz << endl;
  }
#endif

	exit(EXIT_SUCCESS);
}
