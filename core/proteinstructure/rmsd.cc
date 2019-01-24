#include "rmsd.h"

/*
 * This version of the RMSD calculation method assumes that both vectors
 * are the same size and that each position in each array corresponds to
 * the actual alignment of atoms
 */
double calculate_allatom_rmsd(vector<atom>& S, vector<atom>& T)
{
    double rmsd = 0.;
    int N = 0;
    for(size_t atom_index = 0; atom_index < S.size(); atom_index++)
    {
        rmsd += get_squared_distance(S[atom_index].axyz, T[atom_index].axyz);
        N++;
    }
    rmsd /= (double) N;
    rmsd = sqrt(rmsd);
    return rmsd;
}


double calculate_lrms(vector<atom>& S, vector<atom>& T, map<int, int>& E)
{
    double rmsd = 0.;
    int N = 0;
    map<int, int>::iterator iter;
    for(iter = E.begin();iter!=E.end();++iter)
    {
        int q = iter->first;
        int r = iter->second;
        if(S[q].atype != "CA") continue;
        rmsd += get_squared_distance(S[q].axyz, T[r].axyz);
        N++;
    }
    if(N != 0)
    {
        rmsd /= (double) N;
        rmsd = sqrt(rmsd);
    }
    else
    {
        cerr << "No matching CA atoms" << endl;
        exit(EXIT_FAILURE);
    }

    return rmsd;
}

double calculate_irmsd(vector<atom>& BR, vector<atom>& UR, vector<atom>& BL,
    vector<atom>& UL, map<int, int>& MP_R, map<int, int>& MP_L, int k1, int k2)
{
    return align(BR, UR, MP_R, k1, BL, UL, MP_L, k2);
}

// overloaded version that not only calculates fnat fot a target/prediction, but for each of the chains too
// according to the per-chain information given to it
void get_fnat(vector<string>& target, vector<string>& predicted, vector< vector<string> >& target_per_chain, vector< vector<string> >& predicted_per_chain,
			double& fnat, double& fnon_nat, double* fnat_per_chain, double* fnon_nat_per_chain)
{
	// first calculate the overall fnat
	get_fnat(target, predicted, fnat, fnon_nat);
	// then calculate fnat and fnon-nat for each chain
	for(size_t chain_index = 0 ; chain_index < target_per_chain.size(); chain_index++)
	{
		double current_fnat, current_fnon_nat;
		get_fnat(target_per_chain[chain_index], predicted_per_chain[chain_index], current_fnat, current_fnon_nat);
		fnat_per_chain[chain_index] = current_fnat;
		fnon_nat_per_chain[chain_index] = current_fnon_nat;
	}
}

void get_fnat(vector<string>& target, vector<string>& predicted, double& fnat, double& fnon_nat)
{
	double target_contacts = (double) target.size();
	double predicted_contacts = (double) predicted.size();
	double native_contacts = 0.;

	vector<string>::iterator pos;
	string str;
	for(size_t i=0;i<target_contacts;i++)
	{
		pos = find(predicted.begin(), predicted.end(), target[i]);
		if(pos != predicted.end())
		{
			native_contacts += 1.0;
		}
	}

	fnat = target_contacts ? native_contacts / target_contacts : 0;
	fnon_nat = predicted_contacts ? (predicted_contacts - native_contacts) / predicted_contacts : 1;
}

// This function takes a multiple chain target and its prediction and calculates the global c-alpha rmsd
// between them. It directly returns the total rmsd but it also returns the specific rmsd for each chain
// in the output parameter named rmsd_by_chain. The memory for this parameter should have been allocated previously
double calculate_multiple_rmsd(vector< vector<atom> >& target, vector< vector<atom> >& predicted, double* rmsd_by_chain, bool interface_only)
{
	map<pair<int,int>, vector<int> > M;
	// storage for the alignment transformation that will be returned by the alignment
	gsl_matrix *R = gsl_matrix_alloc(3,3); // rotation matrix. 3 by 3 because of the 3 axis values x,y,z
	gsl_vector *t = gsl_vector_alloc(3); // translation vector
	// used to iterate through the pairs of chains
	map<pair<int, int>, vector<int> >::iterator it;
	double squared_error_sum = 0.0;
	int total_ca_contacts = 0;
	int n = get_matching_atoms(predicted, target, M, interface_only);
	double result = -1;
	if(n) {
		result = align(predicted, target, M, n, R, t);
	}

	if(result == -1) { // if -1 is returned no alignment was found because of zero singular values
		for(it=M.begin();it!=M.end();++it) {
			pair<int,int> same_chain_pair = it->first;
			rmsd_by_chain[same_chain_pair.first] = -1;
		}
		return -1;
	}
/*
if(target.size() ==4 && !interface_only) {
	cout << "****** START target\n";
	cerr << "****** START prediction\n";
}*/

	// Iterate through the pairs of chains. There should be a 0-0, 1-1, 2-2, etc...
	for(it=M.begin();it!=M.end();++it)
	{
		double chain_squared_error_sum = 0.0;
		pair<int,int> same_chain_pair = it->first;
		// Get a readable reference to each of the target and predicted chains
		vector<atom> current_predicted_chain = predicted[same_chain_pair.first];
		vector<atom> current_target_chain = target[same_chain_pair.second];
		// This vector has the indices of all the CA atoms that match on both chains
		// We will go through all of these matching atoms and measure the distance
		vector<int> matching_atom_indices = it->second;
		for(size_t atom_index = 0; atom_index < matching_atom_indices.size(); atom_index++)
		{
			int matching_index = matching_atom_indices[atom_index];
			atom current_predicted_atom = current_predicted_chain[matching_index];
			atom current_target_atom = current_target_chain[matching_index];


			// before measuring the distance we need to rotate/translate the predicted
			// points according to the transformation that gives the best aligment
			// described by R and t
			double transformed_predicted_atom[3]; //x,y,z
			// perform the vector X matrix multiplication and add the translation
			// for each of the 3 axis values
			for(size_t axis = 0; axis < 3; axis++)
			{
				transformed_predicted_atom[axis] =
					current_predicted_atom.axyz[0] * gsl_matrix_get(R,axis,0) +
					current_predicted_atom.axyz[1] * gsl_matrix_get(R,axis,1) +
					current_predicted_atom.axyz[2] * gsl_matrix_get(R,axis,2) +
					gsl_vector_get(t,axis);
			}


			chain_squared_error_sum += get_squared_distance(current_target_atom.axyz, transformed_predicted_atom);
		}
/*
if(target.size() ==4 && !interface_only) {
	for(size_t i = 0; i < current_target_chain.size(); i++) {
		cout << current_target_chain[i].tostring();
		atom other = current_predicted_chain[i];
		double transformed_predicted_atom[3];
		for(size_t axis = 0; axis < 3; axis++) {
			transformed_predicted_atom[axis] =
				other.axyz[0] * gsl_matrix_get(R,axis,0) +
				other.axyz[1] * gsl_matrix_get(R,axis,1) +
				other.axyz[2] * gsl_matrix_get(R,axis,2) +
				gsl_vector_get(t,axis);
		}
		other.axyz[0] = transformed_predicted_atom[0];
		other.axyz[1] = transformed_predicted_atom[1];
		other.axyz[2] = transformed_predicted_atom[2];
		cerr << other.tostring();
	}
}*/
		// same_chain_pair contains the indices of the chains. In theory they're both the same
		// The following expression is the RMSD value for the current chain
		rmsd_by_chain[same_chain_pair.first] = sqrt(chain_squared_error_sum / matching_atom_indices.size());
	
		squared_error_sum += chain_squared_error_sum;
		total_ca_contacts += matching_atom_indices.size();
	}
	// release storage from the transformation info
	gsl_matrix_free(R);
	gsl_vector_free(t);
	//return sqrt(squared_error_sum / total_ca_contacts);
	return result;
}
