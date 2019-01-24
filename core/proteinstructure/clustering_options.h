#ifndef _CLUSTERING_OPTIONS_H_
#define _CLUSTERING_OPTIONS_H_

#include <getopt.h>
#include <string>
#include <cstdlib>
#include <vector>
#include <iostream>

using namespace std;


#define INPUT_OPTION "input"
#define INPUT_CHAR 'i'
#define RECEPTOR_CHAR 'R'
#define LIGAND_CHAR 'L'
#define CUTOFF_OPTION "cutoff"
#define CUTOFF_CHAR 'c'
#define USE_LRMSD_OPTION "lrmsd"
#define USE_LRMSD_CHAR 'l'
#define PRINT_SIZES_OPTION "printsizes"
#define PRINT_SIZES_CHAR 'p'


#define SHORT_OPTIONS "R:L:"

static struct option clustering_long_options[] = 
{
	{INPUT_OPTION, required_argument, 0, INPUT_CHAR},
	{CUTOFF_OPTION, required_argument, 0, CUTOFF_CHAR},
	{USE_LRMSD_OPTION, no_argument, 0, USE_LRMSD_CHAR},
	{PRINT_SIZES_OPTION, required_argument, 0, PRINT_SIZES_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the pair-wise LZerD clustering program
class clustering_options
{
	public:
		// Only constructor in this class. It uses the input arguments in order to call getopt_long
		// to help the parsing process. After the instance is created, the methods provided by this class
		// will tell the caller if the process was successful and what is the value of the different configuration
		// parameters used
		clustering_options(int argc, char** argv)
		{
			sizes_file = input_file = receptor_file = ligand_file = "";
			// default value that should be replaced later on
			cutoff = -1;
			// default value for the usage of Ligand RMSD (by default we use all C-alpha atoms)
			use_lrmsd = false;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c; // c represents the character associated with an option
				c = getopt_long(argc, argv, SHORT_OPTIONS, clustering_long_options, &option_index);
				if(c == -1)
				{
					options_left = false;
				}
				else
				{
					// determine what type of option we got
					switch(c)
					{
						case INPUT_CHAR:
							input_file = optarg;
							break;
						case RECEPTOR_CHAR:
							receptor_file = optarg;
							break;
						case LIGAND_CHAR:
							ligand_file = optarg;
							break;
						case CUTOFF_CHAR:
							cutoff = atof(optarg);
							break;
						case USE_LRMSD_CHAR:
							use_lrmsd = true;
							break;
						case PRINT_SIZES_CHAR:
							sizes_file = optarg;
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
					!input_file.empty() &&
					cutoff != -1;
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

		// name of the file with LZerD transformations
		string get_input_file()
		{
			return input_file;
		}

		// If specified, it returns the name of the file where the cluster sizes will be output
		string get_sizes_file()
		{
			return sizes_file;
		}

		// retrieve the flag that determines if we will cluster using C-alpha LRMSD (C-alphas that belong to the ligand
		// within 10 angstroms of at least one receptor atom)
		// or the RMSD of all ligand C-alphas 
		bool get_lrmsd_flag()
		{
			return use_lrmsd;
		}
		
		// retrieve the cutoff (in angstroms) that determines if two elements should be clustered
		// because of their proximity
		double get_cutoff()
		{
			return cutoff;
		}

	private:
		// These variables hold the parameters that come directly from the command line
		string receptor_file;
		string ligand_file;
		string input_file;
		string sizes_file;
		double cutoff;
		bool use_lrmsd;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;

};

#endif
