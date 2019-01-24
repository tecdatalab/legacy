#ifndef _DOCKING_STATS_H_
#define _DOCKING_STATS_H_

using namespace std;

#include <cstdio>

class docking_stats
{
	private:
		double fnat;
		double fnon_nat;
		double irmsd;
		/* Represents global c-alpha rmsd */
		double rmsd;
		/* This is intended to be something like "AB" or "ADF" to describe that the stats correspond to
		 * a subset of the chains being aligned */
		string description;
		/* by default it will be set to zero, when the tostring method is called the idea will be to
		 * output a column that represents the number of subunits that originated this */
		size_t subunits;
		/*
		 * Since these are docking stats, these values should correspond to a measurement from a decoy to
		 * the native. This number establishes which of those decoys it represents
		 */
		int decoy_rank;
	public:
		/*
		 * This is expected to be an immutable object that will get its values from this constructor
		 * and default values can be assigned when necessary
		 */
		docking_stats(double fnat, double fnon_nat, double irmsd, double rmsd, string description="", size_t subunits=0, int decoy_rank=-1) {
			this->fnat = fnat;
			this->fnon_nat = fnon_nat;
			this->irmsd = irmsd;
			this->rmsd = rmsd;
			this->description = description;
			this->subunits = subunits;
			this->decoy_rank = decoy_rank;
		}
		/*
		 * Print a comma separated representation of this stats description
		 * If a true is provided as parameter it will return a string composed of two lines,
		 * the first being a header
		 */
		string to_comma_string(bool print_header=false) {
			string result = "";
			char buf[100];
			if(print_header) {
				if(decoy_rank != - 1) {
					result += "rank,";
				}
				if(subunits != 0) {
					result += "subunits,";
				}
				if(!description.empty()) {
					result += "description,";
				}
				result += "fnat,fnon_nat,irmsd,rmsd\n";
			}
			
			if(decoy_rank != - 1) {
				sprintf(buf, "%d,", decoy_rank);
				result += string(buf);
			}
			if(subunits != 0) {
				sprintf(buf, "%u,", (unsigned int)subunits);
				result += string(buf);
			}
			if(!description.empty()) {
				result += description + ",";
			}
			sprintf(buf, "%f,%f,%f,%f\n", fnat, fnon_nat, irmsd, rmsd);
			result += string(buf);

			return result;
		}
};
#endif

