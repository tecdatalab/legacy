#ifndef _STATSSUBUNITS_H_
#define _STATSSUBUNITS_H_

#include <vector>
#include "atom.h"
#include "docking_stats.h"

using std::string;
using std::vector;

class statssubunits
{
	private:
		vector<vector<atom> > native_complex;
		vector<string> chain_ids;
		vector<vector<size_t> > combinations;
		static vector<vector<size_t> > getindexcombinations(size_t totalsubunits);
		static void docombinations(vector<size_t>* all_indices, vector<size_t>* combinations_template,
						size_t length, size_t level, size_t start, vector<vector<size_t> >* out_combinations);
	public:
		statssubunits(vector<vector<atom> > native_complex, vector<string> chain_ids);
		vector<docking_stats> calculate_stats(vector<vector<atom> > decoy, size_t decoy_rank);
};
#endif

