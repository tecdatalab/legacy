#ifndef _PAIRWISESTATS_OPTIONS_H_
#define _PAIRWISESTATS_OPTIONS_H_

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>

using namespace std;


#define RECEPTOR_CHAR 'R'
#define LIGAND_CHAR 'L'
#define PDB_OPTION "pdb"
#define PDB_CHAR 'p'

#define SHORT_OPTIONS "R:L:"

static struct option pairwisestats_long_options[] = 
{
	{PDB_OPTION, required_argument, 0, PDB_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the program that generates measures to compare
// predictions against the native protein structures (restricted to pair-wise protein-protein predictions)
class pairwisestats_options
{
	public:
		// Only constructor in this class. It uses the input arguments in order to call getopt_long
		// to help the parsing process. After the instance is created, the methods provided by this class
		// will tell the caller if the process was successful and what is the value of the different configuration
		// parameters used
		pairwisestats_options(int argc, char** argv)
		{
			pdb_prediction_file = receptor_file = ligand_file = "";
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c; // c represents the character associated with an option
				c = getopt_long(argc, argv, SHORT_OPTIONS, pairwisestats_long_options, &option_index);
				if(c == -1)
				{
					options_left = false;
				}
				else
				{
					// determine what type of option we got
					switch(c)
					{
						case PDB_CHAR:
							pdb_prediction_file = optarg;
							break;
						case RECEPTOR_CHAR:
							receptor_file = optarg;
							break;
						case LIGAND_CHAR:
							ligand_file = optarg;
							break;
						case '?': //option not recognized (an error is automatically printed)
						default:
							options_left = false;
							parse_error = true;
							break;
					}
				}

			}
		}

		// the parse was successful if the errors flag is false and if all required arguments were supplied
		bool parse_successful()
		{
			return !parse_error && !receptor_file.empty() &&
					!ligand_file.empty() &&
					!pdb_prediction_file.empty();
		}

		// name of the receptor PDB file
		string get_receptor_file()
		{
			return receptor_file;
		}

		// name of the ligand PDB file
		string get_ligand_file()
		{
			return ligand_file;
		}

		// name of the pdb file that will be compared against the native
		string get_pdb_file()
		{
			return pdb_prediction_file;
		}

	private:
		// These variables hold the parameters that come directly from the command line
		string receptor_file;
		string ligand_file;
		string pdb_prediction_file;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;

};

#endif
