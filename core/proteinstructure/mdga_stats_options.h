#ifndef _MDGA_STATS_OPTIONS_H_
#define _MDGA_STATS_OPTIONS_H_

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
#include "md.h"
#include "score_weights_options.h"

using namespace std;


#define INPUT_OPTION "input"
#define INPUT_CHAR 'i'
#define TIMER_OPTION "timer"
#define TIMER_CHAR 't'
#define NO_SCORE_OPTION "noscore"
#define NO_SCORE_CHAR 's'
#define CALPHA_ONLY_OPTION "calphaonly"
#define CALPHA_ONLY_CHAR 'c'
#define BOUND_DIRECTORY_OPTION "boundpath"
#define BOUND_DIRECTORY_CHAR 'b'
#define DIRECTORY_CHAR 'd'
#define DIRECTORY_DEFAULT ""
#define RESULTS_NUMBER_CHAR 'n'


// d is the input directory, where the pdb and prediction files are located
#define SHORT_OPTIONS "d:n:"

static struct option mdga_stats_long_options[] = 
{
	{INPUT_OPTION, required_argument, 0, INPUT_CHAR},
	{WEIGHTS_OPTION, required_argument, 0, WEIGHTS_CHAR},
	{BOUND_DIRECTORY_OPTION, required_argument, 0, BOUND_DIRECTORY_CHAR},
	{TIMER_OPTION, no_argument, 0, TIMER_CHAR},
	{NO_SCORE_OPTION, no_argument, 0, NO_SCORE_CHAR},
	{CALPHA_ONLY_OPTION, no_argument, 0, CALPHA_ONLY_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the multiple docking stats program
class mdga_stats_options : public score_weights_options
{
	public:
		// Only constructor in this class. It uses the input arguments in order to call getopt_long
		// to help the parsing process. After the instance is created, the methods provided by this class
		// will tell the caller if the process was successful and what is the value of the different configuration
		// parameters used
		mdga_stats_options(int argc, char** argv) : score_weights_options()
		{
			input = ""; // this is the only required argument
			directory = DIRECTORY_DEFAULT;
			bound_directory = DIRECTORY_DEFAULT;
			results_number = -1;
			use_timer = false;
			generate_scores = true;
			calpha_only = false;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c; // c represents the character associated with an option
				c = getopt_long(argc, argv, SHORT_OPTIONS, mdga_stats_long_options, &option_index);
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
						case RESULTS_NUMBER_CHAR:
							results_number = atoi(optarg);
							break;
						case WEIGHTS_CHAR:
							this->set_weights(optarg);
                                                        break;
						case TIMER_CHAR:
							use_timer = true;
							break;
						case NO_SCORE_CHAR:
							generate_scores = false;
							break;
						case CALPHA_ONLY_CHAR:
							calpha_only = true;
							break;
						case BOUND_DIRECTORY_CHAR:
							bound_directory = optarg;
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
			return !parse_error && !input.empty() && this->get_weights();
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

		// returns the bound directory if it was supplied or just the original directory if it wasn't
		string get_bound_directory()
		{
			if(bound_directory.empty()) {
				return directory;
			} else {
				return bound_directory;
			}
		}

		// number of results displayed. If -1 return all
		int get_results_number()
		{
			return results_number;
		}
		// true if the total time must be output
		bool get_use_timer()
		{
			return use_timer;
		}
		// true by default, if false only rmsd/fnat values printed
		bool get_generate_scores()
		{
			return generate_scores;
		}
		// false by default
		bool get_calpha_only()
		{
			return calpha_only;
		}
	private:
		// These variables hold the parameters that come directly from the command line further processing is done by the methods
		string directory;
		string input;
		string bound_directory;
		int results_number;
		bool use_timer;
		bool generate_scores;
		bool calpha_only;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;

};

#endif
