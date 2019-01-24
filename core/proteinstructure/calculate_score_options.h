#ifndef _CALCULATE_SCORE_OPTIONS_H_
#define _CALCULATE_SCORE_OPTIONS_H_

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
// Include the command line params configuration for weights
#include "score_weights_options.h"

using namespace std;


#define DECOY_OPTION "decoy"
#define DECOY_CHAR 'd'
#define NATIVE_OPTION "native"
#define NATIVE_CHAR 'n'
#define HYDROGENS_CHAR 'h'

// h is supplied when a pdb with hydrogens added is provided
#define SHORT_OPTIONS "h:"

static struct option calculate_score_long_options[] = 
{
	{DECOY_OPTION, required_argument, 0, DECOY_CHAR},
	{NATIVE_OPTION, required_argument, 0, NATIVE_CHAR},
	{WEIGHTS_OPTION, required_argument, 0, WEIGHTS_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the multiple docking stats program
class calculate_score_options : public score_weights_options
{
	public:
		// Only constructor in this class. It uses the input arguments in order to call getopt_long
		// to help the parsing process. After the instance is created, the methods provided by this class
		// will tell the caller if the process was successful and what is the value of the different configuration
		// parameters used
		calculate_score_options(int argc, char** argv) : score_weights_options()
		{
			decoy = native = hydrogens = "";
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c; // c represents the character associated with an option
				c = getopt_long(argc, argv, SHORT_OPTIONS, calculate_score_long_options, &option_index);
				if(c == -1)
				{
					options_left = false;
				}
				else
				{
					// determine what type of option we got
					switch(c)
					{
						case DECOY_CHAR:
							decoy = optarg;
							break;
						case NATIVE_CHAR:
							native = optarg;
							break;
						case HYDROGENS_CHAR:
							hydrogens = optarg;
							break;
						case WEIGHTS_CHAR:
							this->set_weights(optarg);
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
			return !parse_error && !decoy.empty() && this->get_weights() &&
				((native.empty() && hydrogens.empty()) ||
				 (!native.empty() && !hydrogens.empty()));
		}

		// returns the name of the decoy file
		string get_decoy()
		{
			return decoy;
		}

		// name of the file with hydrogens added
		string get_hydrogens()
		{
			return hydrogens;
		}

		// name of the file with the native conformation
		string get_native()
		{
			return native;
		}

	private:
		// These variables hold the parameters that come directly from the command line further processing is done by the methods
		string decoy;
		string hydrogens;
		string native;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;

};

#endif
