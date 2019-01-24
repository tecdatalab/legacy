#ifndef _MAIN_PDB_ALIGN_OPTIONS_H_
#define _MAIN_PDB_ALIGN_OPTIONS_H_

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>

using namespace std;


#define INPUT_OPTION "input"
#define INPUT_CHAR 'i'
#define OUTPUT_OPTION "output"
#define OUTPUT_CHAR 'o'
#define ALIGN_TO_OPTION "align-to"
#define ALIGN_TO_CHAR 'a'
#define IGNORE_CHAIN_ID_OPTION "ignore-chain-id"
#define IGNORE_CHAIN_ID_CHAR 'c'
#define TRANSFORMED_PDB_OPTION "transformed-pdb"
#define TRANSFORMED_PDB_CHAR 't'

#define SHORT_OPTIONS ""

static struct option main_pdb_align_long_options[] = 
{
	{INPUT_OPTION, required_argument, 0, INPUT_CHAR},
	{OUTPUT_OPTION, required_argument, 0, OUTPUT_CHAR},
	{ALIGN_TO_OPTION, required_argument, 0, ALIGN_TO_CHAR},
	{IGNORE_CHAIN_ID_OPTION, no_argument, 0, IGNORE_CHAIN_ID_CHAR},
	{TRANSFORMED_PDB_OPTION, required_argument, 0, TRANSFORMED_PDB_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the program that
// outputs the atomic coordinates that provide the best alignment between
// the input and the align_to PDB.
class main_pdb_align_options
{
	public:
		main_pdb_align_options(int argc, char** argv)
		{
			input_file = output_file = align_to_file = transformed_file = "";
      ignore_chain_id = false;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c;
				c = getopt_long(argc, argv, SHORT_OPTIONS, main_pdb_align_long_options,
                        &option_index);
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
            case OUTPUT_CHAR:
              output_file = optarg;
              break;
						case ALIGN_TO_CHAR:
							align_to_file = optarg;
							break;
						case IGNORE_CHAIN_ID_CHAR:
							ignore_chain_id = true;
							break;
            case TRANSFORMED_PDB_CHAR:
              transformed_file = optarg;
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
			return !parse_error && !input_file.empty() && !output_file.empty() && !align_to_file.empty();
		}

		string get_input_file()
		{
			return input_file;
		}

    string get_output_file()
    {
      return output_file;
    }

		string get_align_to_file()
		{
			return align_to_file;
		}

    string get_transformed_file()
    {
      return transformed_file;
    }

		bool get_ignore_chain_id()
		{
			return ignore_chain_id;
		}

	private:
		// These variables hold the parameters that come directly from the command line
		string input_file;
    string output_file;
		string align_to_file;
    string transformed_file;
    bool ignore_chain_id;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;
};

#endif
