#include "statssubunits.h"
#include "contact.h"
#include "rmsd.h"

statssubunits::statssubunits(vector<vector<atom> > native_complex, vector<string> chain_ids)
{
	this->native_complex = native_complex;
	this->chain_ids = chain_ids;
	combinations = statssubunits::getindexcombinations(native_complex.size());
}

vector<vector<size_t> > statssubunits::getindexcombinations(size_t totalsubunits)
{
	vector<size_t> combination_template; // this will be used as a buffer to add and remove the indices along the process
	vector<size_t> all_indices; // contains 0,1,2,....,totalsubunits - 1 (it's the possible values for each combination)
	vector<vector<size_t> > all_combinations;
	for(size_t index = 0; index < totalsubunits; index++)
	{
		all_indices.push_back(index);
	}
	statssubunits::docombinations(&all_indices, &combination_template, totalsubunits, 0, 0, &all_combinations);

	return all_combinations;
}

void statssubunits::docombinations(vector<size_t>* all_indices, vector<size_t>* combinations_template,
				size_t length, size_t level, size_t start, vector<vector<size_t> >* out_combinations)
{
	for(size_t i = start; i < length; i++)
	{
		combinations_template->push_back(all_indices->at(i));
		vector<size_t> current_combination = *combinations_template;
		if(current_combination.size() > 1) { // don't include single chains
			out_combinations->push_back(current_combination);
		}

		if(i < length - 1) {
			statssubunits::docombinations(all_indices, combinations_template, length, level + 1, i + 1, out_combinations);
		}
		combinations_template->pop_back();
	}
}

vector<docking_stats> statssubunits::calculate_stats(vector<vector<atom> > decoy, size_t decoy_rank)
{
	/* These 2 vectors will be used to hold the subset of chains analyzed at each step */
	vector<vector<atom> > native_subset;
	vector<vector<atom> > decoy_subset;
	/* Used to set the docking_stats text description */
	string description;
	/* One element for each chain combination in this results collection */
	vector<docking_stats> results;

	// this is required by the method that computes rmsd's but it will not be used
	// it will be initialized to the maximum size possible
	double* rmsd_by_chain = new double[decoy.size()];

	for(size_t combinationindex = 0; combinationindex < combinations.size(); combinationindex++)
	{
		vector<size_t> onecombination = combinations[combinationindex];
		// Add each chain in a combination to both native and decoy subsets
		for(size_t subunitindex = 0; subunitindex < onecombination.size(); subunitindex++)
		{
			size_t subunit = onecombination[subunitindex];
			vector<atom> native_subunit = native_complex[subunit];
			vector<atom> decoy_subunit = decoy[subunit];

			description += chain_ids[subunit];
			if(subunitindex < onecombination.size() - 1) {
				description += "-";
			}

			native_subset.push_back(native_subunit);
			decoy_subset.push_back(decoy_subunit);

		}

		// at this point both subsets should have the appropriate combination and we should compute the statistics
		vector<string> native_contacts;
	        vector< vector<string> > native_contacts_per_chain;
		vector<string> decoy_contacts;
                vector< vector<string> > decoy_contacts_per_chain;
	        identify_contact_residues(native_subset, native_contacts, native_contacts_per_chain);
		identify_contact_residues(decoy_subset, decoy_contacts, decoy_contacts_per_chain);

	        //...as well as marking the interface residues
		get_interface_residues(native_subset);

		// calculate fnat
		double fnat, fnon_nat;
		get_fnat(native_contacts, decoy_contacts, fnat, fnon_nat); //fnat and fnon_nat by ref


		// calculate rmsd
//cout << "******* " << description << endl;
		double rmsd = calculate_multiple_rmsd(native_subset, decoy_subset, rmsd_by_chain, false);
//cout << "RMSD " << rmsd << endl;
		// calculate irmsd
//cout << "+++++++ " << description << endl;
		double irmsd = calculate_multiple_rmsd(native_subset, decoy_subset, rmsd_by_chain, true);
//cout << "iRMSD " << irmsd << endl;

		docking_stats stats(fnat, fnon_nat, irmsd, rmsd, description, native_subset.size(), decoy_rank);
		results.push_back(stats);

		// reset the atoms collections for the next iteration
		native_subset.clear();
		decoy_subset.clear();
		// and the description
		description = "";
	}

	delete [] rmsd_by_chain;

	return results;
}

