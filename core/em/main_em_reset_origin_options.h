#ifndef _MAIN_EM_RESET_ORIGIN_OPTIONS_H_
#define _MAIN_EM_RESET_ORIGIN_OPTIONS_H_

#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>

using namespace std;


#define EM_OPTION "simulated-em"
#define EM_CHAR 'e'
#define PDB_OPTION "base-pdb"
#define PDB_CHAR 'p'
#define CONTOUR_LEVEL_OPTION "contour-level"
#define CONTOUR_LEVEL_CHAR 'c'
#define OUTPUT_OPTION "output"
#define OUTPUT_CHAR 'o'

// In case the user just wants to overwrite the origin using a particular x,y,z
#define SHORT_OPTIONS "x:y:z:"

static struct option main_em_reset_origin_long_options[] = 
{
	{EM_OPTION, required_argument, 0, EM_CHAR},
	{PDB_OPTION, required_argument, 0, PDB_CHAR},
	{CONTOUR_LEVEL_OPTION, required_argument, 0, CONTOUR_LEVEL_CHAR},
	{OUTPUT_OPTION, required_argument, 0, OUTPUT_CHAR},
	{0,0,0,0}
};


// This class parses the command line arguments supplied to the program that
// outputs the atomic coordinates that provide the best alignment between
// the input and the align_to PDB.
class main_em_reset_origin_options
{
	public:
    main_em_reset_origin_options(int argc, char** argv)
		{
      em_file = pdb_file = output_file = "";
      contour_level = 0.0f;
      x_provided = y_provided = z_provided = false;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c;
				c = getopt_long(argc, argv, SHORT_OPTIONS,
                        main_em_reset_origin_long_options, &option_index);
				if(c == -1)
				{
					options_left = false;
				}
				else
				{
					// determine what type of option we got
					switch(c)
					{
						case EM_CHAR:
            	em_file = optarg;
							break;
            case PDB_CHAR:
              pdb_file = optarg;
              break;
						case OUTPUT_CHAR:
							output_file = optarg;
							break;
						case CONTOUR_LEVEL_CHAR:
              contour_level = atof(optarg);
							break;
            case 'x':
              x = atof(optarg);
              x_provided = true;
              break;
            case 'y':
              y = atof(optarg);
              y_provided = true;
              break;
            case 'z':
              z = atof(optarg);
              z_provided = true;
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
			return !parse_error && !em_file.empty() && !output_file.empty()
        && (!pdb_file.empty() || is_direct_origin_provided());
		}

		string get_em_file()
		{
			return em_file;
		}

    string get_pdb_file()
    {
      return pdb_file;
    }

		string get_output_file()
		{
			return output_file;
		}

		float get_contour_level()
		{
			return contour_level;
		}

    float get_x()
    {
      return x;
    }

    float get_y()
    {
      return y;
    }
    
    float get_z()
    {
      return z;
    }

    bool is_direct_origin_provided()
    {
      return x_provided && y_provided && z_provided;
    }

	private:
		// These variables hold the parameters that come directly from the command line
		string em_file;
    string pdb_file;
		string output_file;
    float contour_level;
    float x, y, z;
    bool x_provided, y_provided, z_provided;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;
};

#endif
