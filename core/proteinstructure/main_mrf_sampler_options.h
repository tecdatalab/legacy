#ifndef _MAIN_PDB_ALIGN_OPTIONS_H_
#define _MAIN_PDB_ALIGN_OPTIONS_H_

#include <getopt.h>

#include <cstdlib>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;


#define INPUT_PDBS_CHAR 'i'
#define INPUT_PDBS_DESCRIPTION "n full atom PDB files that are sampled around the centers of the c-alpha\n\t\ttemplate centroids"
#define TEMPLATE_CALPHA_PDB_CHAR 't'
#define TEMPLATE_CALPHA_PDB_DESCRIPTION "c-alpha only PDB template file that represents the best fitting\n\t\tof the units in the complete map. It should contain c-alphas for all subunits in a single file."
#define UNIT_EM_MAPS_CHAR 'e'
#define UNIT_EM_MAPS_DESCRIPTION "n EM map files in MRC format that match each of the subunit PDB files.\n\t\tMust follow the same order as the full atom files."
#define COMPLETE_MAP_OPTION "complete-map"
#define COMPLETE_MAP_CHAR 'c'
#define COMPLETE_MAP_DESCRIPTION "Overall EM map where PDB and EM units are sampled for best possible fit"
#define LABELS_OPTION "labels"
#define LABELS_CHAR 'l'
#define LABELS_DESCRIPTION "Each of the n subunits are assigned a label. This comma-separated list of\n\t\ttext names is matched 1 to 1 in the same order as the full atom PDB files provided."
#define NEIGHBORS_OPTION "neighbors"
#define NEIGHBORS_CHAR 'n'
#define NEIGHBORS_DESCRIPTION "Description of the edges that connect the different units in the graph.\n\t\tThe user must provide these connections."
#define CONTOUR_LEVEL_OPTION "contour-level"
#define CONTOUR_LEVEL_CHAR 'o'
#define CONTOUR_LEVEL_DESCRIPTION "EM map contour level. Usually obtained from the EMDB."

#define SHORT_ROTATION_RANGE_OPTION "short-rotation-range"
#define SHORT_ROTATION_RANGE_CHAR 'a'
#define SHORT_ROTATION_RANGE_DESCRIPTION "Near native decoys are sampled in a range between -short-rotation-range and\n\t\tshort-rotation-range degrees (each axis)"

#define SHORT_ROTATION_STEP_OPTION "short-rotation-step"
#define SHORT_ROTATION_STEP_CHAR 'b'
#define SHORT_ROTATION_STEP_DESCRIPTION "This is the increment in degrees between-short-rotation-range\n\t\tand short-rotation-range"

#define SHORT_TRANSLATION_RANGE_OPTION "short-translation-range"
#define SHORT_TRANSLATION_RANGE_CHAR 'd'
#define SHORT_TRANSLATION_RANGE_DESCRIPTION "Near native decoys are sampled in a range between\n\t\t-short-translation-range and short-translation-range angstroms (each axis)"

#define SHORT_TRANSLATION_STEP_OPTION "short-translation-step"
#define SHORT_TRANSLATION_STEP_CHAR 'f'
#define SHORT_TRANSLATION_STEP_DESCRIPTION "This is the increment in angstroms between -short-translation-range\n\t\tand short-translation-range"

#define LONG_ROTATION_STEP_OPTION "long-rotation-step"
#define LONG_ROTATION_STEP_CHAR 'r'
#define LONG_ROTATION_STEP_DESCRIPTION "Overall rotation space will be sampled using steps of x degrees each,\n\t\tgiven by this parameter."

#define LONG_TRANSLATION_RANGE_OPTION "long-translation-range"
#define LONG_TRANSLATION_RANGE_CHAR 'g'
#define LONG_TRANSLATION_RANGE_DESCRIPTION "Overall translation space is sample between -long-translation-range\n\t\tand long-translation-range angstroms"

#define LONG_TRANSLATION_STEP_OPTION "long-translation-step"
#define LONG_TRANSLATION_STEP_CHAR 'h'
#define LONG_TRANSLATION_STEP_DESCRIPTION "Overall translation space sampling is increment by this amount\n\t\tbetween -long-translation-range and long-translation-range angstroms"

#define OUTPUT_PREFIX_OPTION "output-prefix"
#define OUTPUT_PREFIX_CHAR 'x'
#define OUTPUT_PREFIX_DESCRIPTION "Output filenames will start with this prefix"

#define SWAPPING_ENABLED_OPTION "swapping-enabled"
#define SWAPPING_ENABLED_CHAR 's'
#define SWAPPING_ENABLED_DESCRIPTION "Wrong placement of units are sampled by exchanging them between\n\ttheir centroid locations"

#define PDB_OUTPUT_PREFIX_OPTION "pdb-output-prefix"
#define PDB_OUTPUT_PREFIX_CHAR 'p'
#define PDB_OUTPUT_PREFIX_DESCRIPTION "Triggers the output of PDB files that represent each singleton.\n\tOnly singleton sample files are generated if provided"

#define SHORT_OPTIONS "i:t:e:"

static struct option main_mrf_sampler_long_options[] = 
{
	{COMPLETE_MAP_OPTION, required_argument, 0, COMPLETE_MAP_CHAR},
	{LABELS_OPTION, required_argument, 0, LABELS_CHAR},
	{NEIGHBORS_OPTION, required_argument, 0, NEIGHBORS_CHAR},
	{CONTOUR_LEVEL_OPTION, required_argument, 0, CONTOUR_LEVEL_CHAR},
  {SHORT_ROTATION_RANGE_OPTION, required_argument, 0, SHORT_ROTATION_RANGE_CHAR},
  {SHORT_ROTATION_STEP_OPTION, required_argument, 0, SHORT_ROTATION_STEP_CHAR},
  {SHORT_TRANSLATION_RANGE_OPTION, required_argument, 0,
   SHORT_TRANSLATION_RANGE_CHAR},
  {SHORT_TRANSLATION_STEP_OPTION, required_argument, 0,
   SHORT_TRANSLATION_STEP_CHAR},
  {LONG_ROTATION_STEP_OPTION, required_argument, 0, LONG_ROTATION_STEP_CHAR},
  {LONG_TRANSLATION_RANGE_OPTION, required_argument, 0,
   LONG_TRANSLATION_RANGE_CHAR},
  {LONG_TRANSLATION_STEP_OPTION, required_argument, 0,
   LONG_TRANSLATION_STEP_CHAR},
  {OUTPUT_PREFIX_OPTION, required_argument, 0, OUTPUT_PREFIX_CHAR},
  {PDB_OUTPUT_PREFIX_OPTION, required_argument, 0, PDB_OUTPUT_PREFIX_CHAR},
	{SWAPPING_ENABLED_OPTION, no_argument, 0, SWAPPING_ENABLED_CHAR},
	{0,0,0,0}
};


class main_mrf_sampler_options
{

	public:
		main_mrf_sampler_options(int argc, char** argv)
		{
      complete_map = labels_arg = neighbors_arg = template_pdb = output_prefix = "";
      pdb_output_prefix = "";
      short_rotation_range = short_rotation_step = short_translation_range = 0;
      short_translation_step = long_rotation_step = long_translation_range = 0;
      long_translation_step = 0;
      contour_level = 0;
      swapping_enabled = false;
			// so far we haven't encountered errors
			parse_error = false;

			bool options_left = true;

			while(options_left)
			{
				int option_index, c;
				c = getopt_long(argc, argv, SHORT_OPTIONS, main_mrf_sampler_long_options,
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
						case INPUT_PDBS_CHAR:
							input_pdbs.push_back(optarg);
							break;
            case TEMPLATE_CALPHA_PDB_CHAR:
              template_pdb = optarg;
              break;
						case UNIT_EM_MAPS_CHAR:
							unit_maps.push_back(optarg);
							break;
						case COMPLETE_MAP_CHAR:
              complete_map = optarg;
							break;
            case LABELS_CHAR:
              labels_arg = optarg;
              break;
            case NEIGHBORS_CHAR:
              neighbors_arg = optarg;
              break;
            case CONTOUR_LEVEL_CHAR:
              contour_level = atof(optarg);
              break;
            case SHORT_ROTATION_RANGE_CHAR:
              short_rotation_range = atof(optarg);
              break;
            case SHORT_ROTATION_STEP_CHAR:
              short_rotation_step = atof(optarg);
              break;
            case SHORT_TRANSLATION_RANGE_CHAR:
              short_translation_range = atof(optarg);
              break;
            case SHORT_TRANSLATION_STEP_CHAR:
              short_translation_step = atof(optarg);
              break;
            case LONG_ROTATION_STEP_CHAR:
              long_rotation_step = atof(optarg);
              break;
            case LONG_TRANSLATION_RANGE_CHAR:
              long_translation_range = atof(optarg);
              break;
            case LONG_TRANSLATION_STEP_CHAR:
              long_translation_step = atof(optarg);
              break;
            case OUTPUT_PREFIX_CHAR:
              output_prefix = optarg;
              break;
            case PDB_OUTPUT_PREFIX_CHAR:
              pdb_output_prefix = optarg;
              break;
            case SWAPPING_ENABLED_CHAR:
              swapping_enabled = true;
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
      if(parse_error || complete_map.empty() || labels_arg.empty() ||
         neighbors_arg.empty() || template_pdb.empty() || output_prefix.empty() ||
         long_rotation_step == 0) {
        return false;
      }
      // Check if graph connections belong to label names
      vector<pair<string,string> > neighbors = get_neighbors();
      vector<string> labels = get_labels();
      set<string> labels_set(labels.begin(), labels.end());
      for(vector<pair<string,string> >::iterator it = neighbors.begin();
          it != neighbors.end(); it++) {
        pair<string,string> neighbor_pair = *it;
        if(labels_set.find(neighbor_pair.first) == labels_set.end() ||
           labels_set.find(neighbor_pair.second) == labels_set.end()) {
          return false;
        }
      }
      // Check if the size of the input PDBs and EM maps are the same
      if(input_pdbs.size() != unit_maps.size() ||
         input_pdbs.size() != labels.size()) {
        return false;
      }
      return true;
		}

    string usage()
    {
      return string("Sample usage: mrf_sampler -i A.pdb -i B.pdb -i C.pdb ") +
          "-t ABC-calpha.pdb\n"
          "\t\t-e A.mrf -e B.mrf -e C.mrf --complete-map ABC.mrf "
          "--contour-level 4.07\n"
          "\t\t--labels A,B,C --neighbors A-B,B-C --long-rotation-step 15\n"
          "\t\t--long-translation-range 5 --long-translation-step 1\n"
          "\t\t--short-rotation-range 5 --short-rotation-step 1\n"
          "\t\t--short-translation-range 1 --short-translation-step 0.5\n\n"
          "Parameters:\n\t" +
          string(1, INPUT_PDBS_CHAR) + ": " + INPUT_PDBS_DESCRIPTION + "\n\t" +
          string(1, TEMPLATE_CALPHA_PDB_CHAR) + ": " +
              TEMPLATE_CALPHA_PDB_DESCRIPTION + "\n\t" +
          string(1, UNIT_EM_MAPS_CHAR) + ": " +
              UNIT_EM_MAPS_DESCRIPTION + "\n\t" +
          COMPLETE_MAP_OPTION + ": " + COMPLETE_MAP_DESCRIPTION + "\n\t" +
          LABELS_OPTION + ": " + LABELS_DESCRIPTION + "\n\t" +
          NEIGHBORS_OPTION + ": " + NEIGHBORS_DESCRIPTION + "\n\t" +
          CONTOUR_LEVEL_OPTION + ": " + CONTOUR_LEVEL_DESCRIPTION + "\n\t"
          SHORT_ROTATION_RANGE_OPTION + ": " +
              SHORT_ROTATION_RANGE_DESCRIPTION + "\n\t" +
          SHORT_ROTATION_STEP_OPTION + ": " +
              SHORT_ROTATION_STEP_DESCRIPTION + "\n\t" +
          SHORT_TRANSLATION_RANGE_OPTION + ": " +
              SHORT_TRANSLATION_RANGE_DESCRIPTION + "\n\t" +
          SHORT_TRANSLATION_STEP_OPTION + ": " +
              SHORT_TRANSLATION_STEP_DESCRIPTION + "\n\t" +
          LONG_ROTATION_STEP_OPTION + ": " +
              LONG_ROTATION_STEP_DESCRIPTION + "\n\t" +
          LONG_TRANSLATION_RANGE_OPTION + ": " +
              LONG_TRANSLATION_RANGE_DESCRIPTION + "\n\t" +
          LONG_TRANSLATION_STEP_OPTION + ": " +
              LONG_TRANSLATION_STEP_DESCRIPTION + "\n\t" +
          OUTPUT_PREFIX_OPTION + ": " + OUTPUT_PREFIX_DESCRIPTION + "\n\t" +
          PDB_OUTPUT_PREFIX_OPTION + ": " + PDB_OUTPUT_PREFIX_DESCRIPTION + "\n\t" +
          SWAPPING_ENABLED_OPTION + ": " + SWAPPING_ENABLED_DESCRIPTION + "\n";
    }

		vector<string> get_input_pdbs()
		{
			return input_pdbs;
		}

    string get_template_pdb()
    {
      return template_pdb;
    }

		vector<string> get_unit_maps()
		{
			return unit_maps;
		}

    string get_complete_map()
    {
      return complete_map;
    }

    vector<string> get_labels()
    {
      vector<string> labels;
      string labels_param = labels_arg;
      while(!labels_param.empty())
      {
        size_t separator_pos = labels_param.find_first_of(',');
        string new_label = labels_param.substr(0, separator_pos);
        labels.push_back(new_label);
        labels_param = separator_pos ==  string::npos ?
            "" : labels_param.substr(separator_pos + 1);
      }
      return labels;
    }

    vector<pair<string,string> > get_neighbors()
    {
      vector<pair<string, string> > neighbors;
      string neighbors_param = neighbors_arg;
      while(!neighbors_param.empty())
      {
        size_t separator_pos = neighbors_param.find_first_of(',');
        string neighbors_pair_str = neighbors_param.substr(0, separator_pos);
        size_t hyphen_pos = neighbors_pair_str.find_first_of('-');
        string left = neighbors_pair_str.substr(0, hyphen_pos);
        string right = neighbors_pair_str.substr(hyphen_pos + 1);
        neighbors.push_back(make_pair(left,right));

        neighbors_param = separator_pos ==  string::npos ?
            "" : neighbors_param.substr(separator_pos + 1);
      }
      return neighbors;
    }

    double get_contour_level()
    {
      return contour_level;
    }

    double get_short_rotation_range()
    {
      return short_rotation_range;
    }

    double get_short_rotation_step()
    {
      return short_rotation_step;
    }

    double get_short_translation_range()
    {
      return short_translation_range;
    }

    double get_short_translation_step()
    {
      return short_translation_step;
    }

    double get_long_rotation_step()
    {
      return long_rotation_step;
    }

    double get_long_translation_range()
    {
      return long_translation_range;
    }

    double get_long_translation_step()
    {
      return long_translation_step;
    }

    string get_output_prefix()
    {
      return output_prefix;
    }

    string get_pdb_output_prefix()
    {
      return pdb_output_prefix;
    }

    bool get_swapping_enabled()
    {
      return swapping_enabled;
    }

	private:
    vector<string> input_pdbs;
    string template_pdb;
    vector<string> unit_maps;
    string complete_map;
    string labels_arg;
    string neighbors_arg;
    string output_prefix;
    string pdb_output_prefix;
    double short_rotation_range, short_rotation_step, short_translation_range;
    double short_translation_step, long_rotation_step, long_translation_range;
    double long_translation_step;
    double contour_level;
    bool swapping_enabled;
		// set to false in case that unrecognized options or errors are found by getopt_long
		bool parse_error;
};

#endif
