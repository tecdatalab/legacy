#ifndef _MD_H_
#define _MD_H_

#include "pdb.h"
#include "transformations.h"
#include "zdock_transformations.h"
#include "lzerd_transformations.h"

#define CHAIN_FILE_SEPARATOR "-"
#define PDB_EXTENSION ".pdb"
#define PDB_WITH_HYDROGEN_EXTENSION ".pdb.h"
#define PREDICTION_EXTENSION ".out"

// Standard string used to represent decoy programs and scoring types
#define SCORETYPE_PAIRWISE "pairwise"
#define SCORETYPE_PHYSICS "physics"
#define DECOYPROGRAM_LZERD "LZerD"
#define DECOYPROGRAM_ZDOCK "ZDOCK"

/*
 * Defines the types of scoring handled by the multiple docking program, namely, using the
 * averaged pairwise score or the physics based approach
 */
enum score_type_t {score_type_pairwise, score_type_physics};

/*
 * The multiple docking program takes two types of input: zdock and LZerD decoys
 */
enum decoy_program_t {decoy_program_lzerd, decoy_program_zdock};

/*
 * Since we should have generated one file per chain, we assume that there are
 * several files named following the <complex_name>-<chain_id>.pdb pattern
 * This function calls the read_protein function in pdb.cc using each
 * of the filenames that should exist.
 * use_hydrogens indicates if we should load the .pdb.h files that have hydrogens added to them
 * It returns a vector of all the pdb instances loaded.
 */
vector<pdb> load_pdbs(string main_complex_name, vector<string> chain_ids, string directory, bool use_hydrogens, bool calpha_only=false);

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
vector< vector<transformations*> > load_predictions(vector<string> chain_ids, string directory, decoy_program_t decoy_program);

#endif
