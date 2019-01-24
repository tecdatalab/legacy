#include "md.h"

using namespace std;

#include <iostream>

/*
 * Since we should have generated one file per chain, we assume that there are
 * several files named following the <complex_name>-<chain_id>.pdb pattern
 * This function calls the read_protein function in pdb.cc using each
 * of the filenames that should exist.
 * use_hydrogens indicates if we should load the .pdb.h files that have hydrogens added to them
 * It returns a vector of all the pdb instances loaded.
 */
vector<pdb> load_pdbs(string main_complex_name, vector<string> chain_ids, string directory, bool use_hydrogens, bool calpha_only)
{
	// create a pdb container for each chain
	vector<pdb> all_proteins(chain_ids.size(), pdb());
	for(size_t chain = 0; chain < chain_ids.size(); chain++)
	{
		string extension = use_hydrogens ? PDB_WITH_HYDROGEN_EXTENSION : PDB_EXTENSION;
		string chain_filename = directory + chain_ids[chain] + CHAIN_FILE_SEPARATOR +
									main_complex_name + extension;
		read_protein(chain_filename.c_str(), all_proteins[chain], calpha_only);
	}
	return all_proteins;
}

/*
 * Assuming that there are several <chain_id>-<chain_id>.out files, that
 * contain ZDock predictions, this function loads the information into zdata
 * instances, using read_zdock_data function in zdock.cc
 * The order in which chain ids are supplied is important, because it assumes
 * a specific sequence in which predictions have been generated. For instance,
 * if we have chains A, B, C and D, then we would expect files for A-B, A-C,
 * A-D, B-C, B-D and C-D.
 * It is necessary to specify if the decoys where generated with ZDOCK or LZerD
 */
vector< vector<transformations*> > load_predictions(vector<string> chain_ids, string directory, decoy_program_t decoy_program)
{
	// initialize a matrix of predictions, setting the values to NULL
	// since we are not having predictions for all combinations
	// we will invoke the "invert" method to create opposite transformations
	// with respect to the ones that we did calculate using the decoy programs
	vector< vector<transformations*> > all_predictions(chain_ids.size(),
						vector<transformations*>(chain_ids.size(), NULL));
	
	for(size_t receptor_index = 0; receptor_index < chain_ids.size(); receptor_index ++)
	{
		// for the current receptor, go from receptor_index + 1 up to N
		for(size_t ligand_index = receptor_index + 1; ligand_index < chain_ids.size();
				ligand_index++)
		{
			string prediction_filename = directory + chain_ids[receptor_index] +
								CHAIN_FILE_SEPARATOR +
								chain_ids[ligand_index] + PREDICTION_EXTENSION;
			/* this will hold the ligand/receptor transformation, obtained by invoking the invert method */
			transformations* inverted;
			if(decoy_program == decoy_program_zdock) // then create zdock_transformations instances
			{
				all_predictions[receptor_index][ligand_index] = new zdock_transformations(prediction_filename);
				inverted = new zdock_transformations(prediction_filename);
			}
			else // lzerd_transformations instance
			{
				all_predictions[receptor_index][ligand_index] = new lzerd_transformations(prediction_filename);
				inverted = new lzerd_transformations(prediction_filename);
			}
			/* apply the inversion and store it */
			inverted->invert();

			all_predictions[ligand_index][receptor_index] = inverted;
		}
	}
	return all_predictions;
}
