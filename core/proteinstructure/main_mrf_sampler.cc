#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "main_mrf_sampler_options.h"
#include "mrf_sampler.h"

int main(int argc, char** argv)
{
	main_mrf_sampler_options options(argc,argv);
	if(!options.parse_successful())
	{
    std::cerr << options.usage();
		exit(EXIT_FAILURE);
	}

  vector<string> input_pdb_files = options.get_input_pdbs();
  string calpha_template_file = options.get_template_pdb();
  vector<string> map_files = options.get_unit_maps();
  string complete_map_file = options.get_complete_map();
  vector<string> labels = options.get_labels();
  vector<pair<string,string> > neighbors = options.get_neighbors();
  double contour_level = options.get_contour_level();
  double short_rotation_range = options.get_short_rotation_range();
  double short_rotation_step = options.get_short_rotation_step();
  double short_translation_range = options.get_short_translation_range();
  double short_translation_step = options.get_short_translation_step();
  double long_rotation_step = options.get_long_rotation_step();
  double long_translation_range = options.get_long_translation_range();
  double long_translation_step = options.get_long_translation_step();
  string output_prefix = options.get_output_prefix();
  string pdb_output_prefix = options.get_pdb_output_prefix();
  double swapping_enabled = options.get_swapping_enabled();

  std::cout << "Sampling parameters:\n"
            << "\tPrefix: " << output_prefix << endl
            << "\tShort rotation plus-minus: " << short_rotation_range << endl
            << "\tShort rotation step: " << short_rotation_step << endl
            << "\tShort translation plus-minus: " << short_translation_range
            << endl
            << "\tShort translation step: " << short_translation_step << endl
            << "\tLong rotation step: " << long_rotation_step << endl
            << "\tLong translation plus-minus: " << long_translation_range
            << endl
            << "\tLong translation step: " << long_translation_step << endl;

  mrf_sampler sampler(input_pdb_files, calpha_template_file, map_files,
                      complete_map_file, contour_level, labels, neighbors,
                      short_rotation_range, short_rotation_step,
                      short_translation_range, short_translation_step,
                      long_rotation_step, long_translation_range,
                      long_translation_step, output_prefix, pdb_output_prefix);


  sampler.generate_feature_files(swapping_enabled);

	exit(EXIT_SUCCESS);
}
