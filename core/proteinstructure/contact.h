#ifndef _CONTACT_H_
#define _CONTACT_H_

#include "pdb.h"
#include "align.h"
#include <set>
#include <algorithm>
#include <ANN/ANN.h>

void get_interface_residues(vector<atom>&, vector<atom>&);
void get_interface_residues(vector< vector<atom> >& proteins);
void set_ligand_CA(map<int, int>&, vector<string>&, vector<atom>&, vector<atom>&, int&);
void check_ligand_CA(map<int, int>&);
void identify_contact_residues(vector<atom>&, vector<atom>&, ANNkd_tree*, vector<string>&);
void get_residue_list(vector<atom>&, vector<string>&);

void check_CA(map<int, int>&, int);
void set_CA(map<int, int>&, vector<atom>&, vector<atom>&, int&);
void identify_contact_residues(vector<atom>&, vector<atom>&, ANNkd_tree*, vector<string>&);
void identify_contact_residues(vector< vector<atom> > proteins, vector<string>& MP, vector< vector<string> >& contacts_per_chain);
int get_matching_atoms(vector<vector<atom> >& A, vector<vector<atom> >& B, map<pair<int,int>, vector<int> >& M,bool interface_only=false);

#endif
