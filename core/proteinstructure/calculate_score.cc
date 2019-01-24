#include <iostream>
#include <vector>
#include "pdb.h"
#include "atom.h"
#include "scoring.h"
#include "calculate_score_options.h"
#include "align.h"
#include "contact.h"
#include "md.h"
#include "mdgraph.h"
#include "rmsd.h"
#include "soroban_score.h"
#include "transformations.h"
#include "utils.h"

using std::cerr;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char** argv)
{
	calculate_score_options options(argc,argv);
	if(!options.parse_successful())
	{
		cerr << "Usage: score --decoy <pdb file> [-h <pdb file + hydrogens> --native <native pdb file>] " << WEIGHT_PARAM_DESCRIPTION << endl;
		exit(EXIT_FAILURE);
	}

	// Generate RMSD's only if we are provided with a native file
	bool generate_distances = false;


	string decoy_file = options.get_decoy();
	string native_file = options.get_native();
	string hydrogens_file = options.get_hydrogens();

	double* weights = options.get_weights();

	// each possible protein complex. Since some of them are optional they are set to empty first
	vector<vector<atom> > decoy = read_chains(decoy_file);
	vector<vector<atom> > native, hydrogens;
	if(!native_file.empty())
	{
		native = read_chains(native_file);
		generate_distances = true;
	}
	if(!hydrogens_file.empty())
	{
		hydrogens = read_chains(hydrogens_file);
	}

	// print RMSD statistics if a native conformation was provided
	if(generate_distances)
	{
		// generate the headers first using the chain names
		// assume that each chain has at least one atom and that it has the chain set (it should)
		string fnatheaders, rmsdheaders, irmsdheaders;
		for(size_t chain_index = 0; chain_index < decoy.size(); chain_index++)
		{
			string current_chain = decoy[chain_index][0].chain;

			fnatheaders += ",fnat " + current_chain + ",fnon-nat " + current_chain;
			rmsdheaders += ",RMSD " + current_chain;
			irmsdheaders += ",iRMSD " + current_chain;
		}
		// print all headers
		cout << "fnat,fnon-nat" << fnatheaders << ",Avg fnat,Avg fnon-nat,RMSD" << rmsdheaders << ",Avg RMSD,iRMSD" <<
			irmsdheaders << ",Avg iRMSD,Clashes,Score" <<
			",vdw,vdw_attr,vdw_rep,elec,elec_sr_attr,elec_lr_attr,elec_sr_rep,elec_lr_rep,hbp_ss,solv,sasa,acp" << endl;

		// get the contact information
		vector<string> original_contacts;
		vector< vector<string> > original_contacts_per_chain;
		identify_contact_residues(native, original_contacts, original_contacts_per_chain);
		//...as well as marking the interface residues
		get_interface_residues(native);

		//This array will hold the rmsd values for each chain
		double* rmsd_by_chain = new double[decoy.size()];

		// average variables
		double avg_fnat, avg_fnon_nat, avg_rmsd, avg_irmsd;
		// variable used for fnat storage, they are reset on each iteration
		double fnat, fnon_nat;
		double* fnat_per_chain = new double[decoy.size()];
		double* fnon_nat_per_chain = new double[decoy.size()];

		// set the averages to zero
		avg_fnat = avg_fnon_nat = avg_rmsd = avg_irmsd = 0;

		vector<string> individual_contacts;
		vector< vector<string> > individual_contacts_per_chain;

		// calculate fnat
		identify_contact_residues(decoy, individual_contacts, individual_contacts_per_chain);
		get_fnat(original_contacts, individual_contacts, original_contacts_per_chain, individual_contacts_per_chain,
				fnat, fnon_nat, fnat_per_chain, fnon_nat_per_chain);
		//...and then the RMSD
		double rmsd = calculate_multiple_rmsd(native, decoy, rmsd_by_chain, false);

		cout << fnat << "," << fnon_nat;
		// print the fnat/fnon-nat for each chain
		for(size_t fnat_index = 0; fnat_index < native.size(); fnat_index++)
		{
			// update averages
			avg_fnat += fnat_per_chain[fnat_index];
			avg_fnon_nat += fnon_nat_per_chain[fnat_index];

			cout << "," << fnat_per_chain[fnat_index] << "," << fnon_nat_per_chain[fnat_index];
		}

		// divide by the number of chains
		avg_fnat /= native.size();
		avg_fnon_nat /= native.size();
		cout << "," << avg_fnat << "," << avg_fnon_nat;

		cout << "," << rmsd;
		// print the rmsd for each chain
		for(size_t rmsd_index = 0; rmsd_index < native.size(); rmsd_index++)
		{
			// update average
			avg_rmsd += rmsd_by_chain[rmsd_index];
			cout << "," << rmsd_by_chain[rmsd_index];
		}

		// divide by the number of chains
		avg_rmsd /= native.size();
		cout << "," << avg_rmsd;

		// as well as the irmsd
		double irmsd = calculate_multiple_rmsd(native, decoy, rmsd_by_chain, true);
		cout << "," << irmsd;
		for(size_t rmsd_index = 0; rmsd_index < native.size(); rmsd_index++)
		{
			// update average
			avg_irmsd += rmsd_by_chain[rmsd_index];
			cout << "," << rmsd_by_chain[rmsd_index];
		}
		// divide by the number of chains
		avg_irmsd /= native.size();

		// print a final comma because the soroban score is going to be output next
		cout << "," << avg_irmsd << "," << get_clashing_atoms(decoy).size() << ",";

		delete[] rmsd_by_chain;
		delete[] fnat_per_chain;
		delete[] fnon_nat_per_chain;
	}
	else
	{
		// no rmsd values so just print the score headers
		cout<< "Score,vdw,vdw_attr,vdw_rep,elec,elec_sr_attr,elec_lr_attr,elec_sr_rep,elec_lr_rep,hbp_ss,solv,sasa,acp" << endl;
	}

	soroban_score score = hydrogens.empty() ? compute_energy(decoy, weights) : compute_energy(hydrogens, weights);
	cout << score.calculate_weighted_score() << ",";
	score.print(cout);
	cout << endl;

	exit(EXIT_SUCCESS);
}
