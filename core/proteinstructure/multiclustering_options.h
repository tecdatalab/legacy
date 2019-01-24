#ifndef _MULTICLUSTERING_OPTIONS_H_
#define _MULTICLUSTERING_OPTIONS_H_

#include <getopt.h>
#include <string>
#include <cstdlib>
#include <vector>
#include <iostream>

using namespace std;

#define CUTOFF_OPTION "cutoff"
#define CUTOFF_CHAR 'c'
#define PRINT_SIZES_OPTION "printsizes"
#define PRINT_SIZES_CHAR 'p'
#define INPUT_OPTION "input"
#define INPUT_CHAR 'i'
#define DIRECTORY_CHAR 'd'
#define DIRECTORY_DEFAULT ""
#define OUTPUT_CHAR 'o'


// d is the input directory, where the pdb and prediction files are located
#define MULTICLUSTERING_SHORT_OPTIONS "d:o:"

static struct option multiclustering_long_options[] = 
{
	{INPUT_OPTION, required_argument, 0, INPUT_CHAR},
	{CUTOFF_OPTION, required_argument, 0, CUTOFF_CHAR},
	{PRINT_SIZES_OPTION, required_argument, 0, PRINT_SIZES_CHAR},
	{0,0,0,0}
};

// This class parses the command line arguments supplied to the multiple docking clustering program
class multiclustering_options
{
	public:
		// Only constructor in this class. It uses the input arguments in order to call getopt_long
		// to help the parsing process. After the instance is created, the methods provided by this class
		// will tell the caller if the process was successful and what is the value of the different configuration
		// parameters used
		multiclustering_options(int argc, char** argv)
		{
			input = ""; // this is the only required argument
			directory = DIRECTORY_DEFAULT;
			
			sizes_file = "";
			// default value that should be replaced later on
			cutoff = -1;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c; // c represents the character associated with an option
				c = getopt_long(argc, argv, MULTICLUSTERING_SHORT_OPTIONS, multiclustering_long_options, &option_index);
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
							input = optarg;
							break;
						case DIRECTORY_CHAR:
							directory = optarg;
							break;
						case CUTOFF_CHAR:
							cutoff = atof(optarg);
							break;
						case PRINT_SIZES_CHAR:
							sizes_file = optarg;
							break;
						case OUTPUT_CHAR:
							output_file = optarg;
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
			return !parse_error && !input.empty() && !output_file.empty() && cutoff != -1;
		}

		// returns the input ga file provided
		string get_input_file()
		{
			return input;
		}

		// returns the name of the input directory provided or its default value, if it wasn't provided
		string get_input_directory()
		{
			return directory;
		}

		/* If specified, it returns the name of the file where the cluster sizes will be output */
		string get_sizes_file()
		{
			return sizes_file;
		}

		/*
		 * Retrieve the cutoff (in angstroms) that determines if two elements should be clustered
		 * because of their proximity
		 */
		double get_cutoff()
		{
			return cutoff;
		}

		/*
		 * Return the name of the output file that will be used
		 */
		string get_output_file()
		{
			return output_file;
		}
	
	private:
		// These variables hold the parameters that come directly from the command line further processing is done by the methods
		string directory;
		string input;
		string sizes_file;
		string output_file;
		double cutoff;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;
};

#endif
