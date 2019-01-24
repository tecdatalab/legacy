#ifndef _MAIN_EM_CORRELATION_OPTIONS_H_
#define _MAIN_EM_CORRELATION_OPTIONS_H_

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#define REFERENCE_OPTION "reference"
#define REFERENCE_CHAR 'r'
#define REFERENCE_DESCRIPTION "File that holds the big MRC map that units will be fit into"
#define FITTED_OPTION "fitted"
#define FITTED_CHAR 'f'
#define FITTED_DESCRIPTION "File that holds the MRC map with the subunit"
#define CONTOUR_OPTION "contour"
#define CONTOUR_CHAR 'c'
#define CONTOUR_DESCRIPTION "Density value threshold to create the surfaces"
#define CORRELATION_OPTION "correlation"
#define CORRELATION_CHAR 'n'
#define CORRELATION_DESCRIPTION "Only predictions with CC over this value are returned"
#define OVERLAP_OPTION "overlap"
#define OVERLAP_CHAR 'o'
#define OVERLAP_DESCRIPTION "Only predictions with overlap higher than this are returned"
#define STEP_OPTION "step"
#define STEP_CHAR 's'
#define STEP_DESCRIPTION "Number of voxels shifted every time a new pose is evaluated"
#define ROTATION_STEP_OPTION "rotation"
#define ROTATION_STEP_CHAR 't'
#define ROTATION_STEP_DESCRIPTION "Rotation step size (in degrees)"
#define OUTPUT_OPTION "output"
#define OUTPUT_CHAR 'u'
#define OUTPUT_DESCRIPTION "File where results will be stored"

// No short options in this binary
#define SHORT_OPTIONS ""

static struct option main_em_correlation_long_options[] = 
{
	{REFERENCE_OPTION, required_argument, 0, REFERENCE_CHAR},
	{FITTED_OPTION, required_argument, 0, FITTED_CHAR},
	{CONTOUR_OPTION, required_argument, 0, CONTOUR_CHAR},
	{CORRELATION_OPTION, required_argument, 0, CORRELATION_CHAR},
	{OVERLAP_OPTION, required_argument, 0, OVERLAP_CHAR},
	{STEP_OPTION, required_argument, 0, STEP_CHAR},
	{ROTATION_STEP_OPTION, required_argument, 0, ROTATION_STEP_CHAR},
	{OUTPUT_OPTION, required_argument, 0, OUTPUT_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the the em
// correlation calculation program
class main_em_correlation_options
{
    private:
        // These variables hold the parameters that come directly from the command line further processing is done by the methods
        string reference;
        string fitted;
        float contour;
		float correlation;
		float overlap;
        int step;
		float rotation_step;
		string output_filename;
        // set to false in case that unrecognized options or errors are found by getopt_long
        bool parse_error;

	public:
		// Only constructor in this class. It uses the input arguments in order to call getopt_long
		// to help the parsing process. After the instance is created, the methods provided by this class
		// will tell the caller if the process was successful and what is the value of the different configuration
		// parameters used
		main_em_correlation_options(int argc, char** argv)
		{
			reference = fitted = output_filename = "";
			// default values, it's OK if they're not provided
            contour = correlation = overlap = 0; 
            step = 1;
			rotation_step = 90;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c; // c represents the character associated with an option
				c = getopt_long(argc, argv, SHORT_OPTIONS, main_em_correlation_long_options, &option_index);
				if(c == -1)
				{
					options_left = false;
				}
				else
				{
					// determine what type of option we got
					switch(c)
					{
						case REFERENCE_CHAR:
                            reference = optarg;
							break;
						case FITTED_CHAR:
							fitted = optarg;
							break;
						case CONTOUR_CHAR:
							contour = atof(optarg);
							break;
						case CORRELATION_CHAR:
							correlation = atof(optarg);
							break;
						case OVERLAP_CHAR:
							overlap = atof(optarg);
							break;
						case STEP_CHAR:
							step = atoi(optarg);
							break;
						case ROTATION_STEP_CHAR:
							rotation_step = atof(optarg);
							break;
						case OUTPUT_CHAR:
							output_filename = optarg;
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

		// the parse was successful if the errors flag is false and if all required
        // arguments were supplied
		bool parse_successful()
		{
			return !parse_error && !reference.empty() && !fitted.empty()
					&& !output_filename.empty();
		}
    
        string usage()
        {
            return
                string("Usage: em_correlation --reference mrc_file --fitted mrc_file ") +
                "[--contour density_value] [--correlation cc_threshold] " +
      					"[--overlap overlap_fraction] [--step voxel_step_size] " + 
			      		"[--rotation rotation_step_size] --output filename\n\n" +
                "Parameters:\n\t" +
                REFERENCE_OPTION + ": " + REFERENCE_DESCRIPTION + "\n\t" +
                FITTED_OPTION + ": " + FITTED_DESCRIPTION + "\n\t" +
                CONTOUR_OPTION + ": " + CONTOUR_DESCRIPTION + "\n\t" +
                CORRELATION_OPTION + ": " + CORRELATION_DESCRIPTION + "\n\t" +
                OVERLAP_OPTION + ": " + OVERLAP_DESCRIPTION + "\n\t" +
                STEP_OPTION + ": " + STEP_DESCRIPTION + "\n\t" +
					      ROTATION_STEP_OPTION + ": " + ROTATION_STEP_DESCRIPTION + "\n\t" +
      					OUTPUT_OPTION + ": " + OUTPUT_DESCRIPTION +
			      		"\n";
        }

		// returns the name of the reference map file
		string get_reference_map()
		{
			return reference;
		}

		// returns the name of the fitted map file
		string get_fitted_map()
		{
			return fitted;
		}

		// contour value that defines the appropriate map surface
		float get_contour_value()
		{
			return contour;
		}

		// Minimum CC satisfied to include the transformation
		float get_correlation_threshold()
		{
			return correlation;
		}

		// Minimum overlap required to be included in the results
		float get_overlap_threshold()
		{
			return overlap;
		}
    
        // Number of voxels moved after
        int get_voxel_step_size()
        {
            return step;
        }
		// Angle, in degrees, that serves as rotation step
		float get_rotation_step_size()
		{
			return rotation_step;
		}
		// Path to the output file created
		string get_output_filename()
		{
			return output_filename;
		}
};

#endif
