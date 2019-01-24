#ifndef _MONTECARLO_OPTIONS_H_
#define _MONTECARLO_OPTIONS_H_

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
#include "md.h"

using namespace std;


#define INPUT_OPTION "input"
#define INPUT_CHAR 'i'
#define TRIALS_OPTION "trials"
#define TRIALS_CHAR 't'
#define PREFIX_OPTION "prefix"
#define PREFIX_CHAR 'p'
#define PREFIX_DEFAULT "refined"
#define DIRECTORY_CHAR 'd'
#define DIRECTORY_DEFAULT ""
#define RESULTS_NUMBER_CHAR 'n'


// d is the input directory, where the pdb and prediction files are located
#define SHORT_OPTIONS "d:n:"

static struct option montecarlo_long_options[] = 
{
	{INPUT_OPTION, required_argument, 0, INPUT_CHAR},
	{TRIALS_OPTION, required_argument, 0, TRIALS_CHAR},
	{PREFIX_OPTION, required_argument, 0, PREFIX_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the montecarlo refinement program
class montecarlo_options
{
	public:
		// Only constructor in this class. It uses the input arguments in order to call getopt_long
		// to help the parsing process. After the instance is created, the methods provided by this class
		// will tell the caller if the process was successful and what is the value of the different configuration
		// parameters used
		montecarlo_options(int argc, char** argv)
		{
			input = ""; // this is the only required argument
			directory = DIRECTORY_DEFAULT;
			prefix = PREFIX_DEFAULT;
			results_number = -1;
			trials = -1;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c; // c represents the character associated with an option
				c = getopt_long(argc, argv, SHORT_OPTIONS, montecarlo_long_options, &option_index);
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
						case PREFIX_CHAR:
							prefix = optarg;
							break;
						case RESULTS_NUMBER_CHAR:
							results_number = atoi(optarg);
							break;
						case TRIALS_CHAR:
							trials = atoi(optarg);
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
			return !parse_error && !input.empty() && trials != -1;
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
		// pdb prefix used to create the actual pdb files
		string get_prefix()
		{
			return prefix;
		}

		// number of results displayed. If -1 return all
		int get_results_number()
		{
			return results_number;
		}
		// number if trials performed for each montecarlo optimization
		int get_trials()
		{
			return trials;
		}
		
	private:
		// These variables hold the parameters that come directly from the command line further processing is done by the methods
		string directory;
		string input;
		string prefix;
		int results_number;
		int trials;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;

};

#endif
