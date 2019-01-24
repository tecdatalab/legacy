#ifndef _RMSD_H_
#define _RMSD_H_

#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_matrix.h>
#include "pdb.h"
#include "atom.h"
#include "align.h"
#include "contact.h"

// Indicates that, when calculating the ligand rmsd between a receptor/ligand pair
// only atoms (in the ligand) closer than 10 angstroms to any receptor atom should
// be considered
#define LIGAND_RMSD_THRESHOLD 10

/*
 * This version of the RMSD calculation method assumes that both vectors
 * are the same size and that each position in each array corresponds to
 * the actual alignment of atoms
 */
double calculate_allatom_rmsd(vector<atom>& S, vector<atom>& T);
double calculate_lrms(vector<atom>&, vector<atom>&, map<int, int>&);
double calculate_irmsd(vector<atom>&, vector<atom>&, vector<atom>&,
    vector<atom>&, map<int, int>&, map<int, int>&, int, int);
void get_fnat(vector<string>& target, vector<string>& predicted, double& fnat, double& fnon_nat);
// overloaded version that not only calculates fnat fot a target/prediction, but for each of the chains too
// according to the per-chain information given to it
void get_fnat(vector<string>& target, vector<string>& predicted, vector< vector<string> >& target_per_chain, vector< vector<string> >& predicted_per_chain,
                         double& fnat, double& fnon_nat, double* fnat_per_chain, double* fnon_nat_per_chain);
// This function takes a multiple chain target and its prediction and calculates the global c-alpha rmsd
// between them. It directly returns the total rmsd but it also returns the specific rmsd for each chain
// in the output parameter named rmsd_by_chain. The memory for this parameter should have been allocated previously
double calculate_multiple_rmsd(vector< vector<atom> >& target, vector< vector<atom> >& predicted, double* rmsd_by_chain, bool interface_only);
#endif
