#include "density_lattice.h"
#include "main_em_correlation_options.h"
#include "mrc.h"

#include <cstdlib>
#include <iostream>

#ifdef WITH_MPI

// only include if the mpi version is being compiled
#include "mpi.h"

#endif


/*
 * Loads MRC files and computes the correlation coefficients
 * between them using different alignments
 */
int main(int argc, char** argv)
{
#ifdef WITH_MPI
    /*** Perform the MPI Initialization that applies to all ranks ***/
    MPI_Init(&argc,&argv);
#endif

	main_em_correlation_options options(argc,argv);
  if(!options.parse_successful()) {
  	cout << options.usage();
    return 0;
  }
  mrc_reader* complete_reader = new mrc_reader(options.get_reference_map());
  mrc_reader* section_reader = new mrc_reader(options.get_fitted_map());
  float density_threshold = options.get_contour_value();
	float correlation_threshold = options.get_correlation_threshold();
	float overlap_threshold = options.get_overlap_threshold();
  int step_size = options.get_voxel_step_size();
	double rotation_step_size = options.get_rotation_step_size();
	string output_filename = options.get_output_filename();
  density_lattice* complete_densities = complete_reader->get_densities();
  density_lattice* section_densities = section_reader->get_densities();

	correlation_results results =
		section_densities->get_significant_correlations(complete_densities,
			step_size, rotation_step_size, density_threshold,
			correlation_threshold, overlap_threshold);
#ifdef WITH_MPI
  // To avoid write conflicts, if we're using MPI write
  // each rank's results to a separate file.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	char rank_number_str[4];
	sprintf(rank_number_str, "%03i", rank);
	output_filename = output_filename + string("-") + rank_number_str ;
#endif
	results.write(output_filename);
/*
	cout << "#X-trans Y-trans Z-trans X-rot Y-rot Z-rot CC Overlap" << endl;
	for(size_t i=0; i< results.size(); i++) {
		correlation_result r = results[i];
	    cout << r.x << " " << r.y << " " << r.z  << " " << r.gamma_x << " "
					 << r.beta_y << " " << r.alpha_z << " "
                     << r.correlation << " " << r.overlap << endl;
	}*/
/*    for (unsigned long x = 0; x <= 100; x+=step_size) {
        for (unsigned long y = 0; y <= 100; y+=step_size) {
            for (unsigned long z = 0; z <= 100; z+=step_size) {
                unsigned long overlap = 0;
                double cc =
                    complete_densities->calculate_correlation(section_densities,
                                                             x, y, z,
                                                             density_threshold,
                                                             &overlap);
                if (cc != 0.0) {
                    cout << "X:" << x << " Y:" << y << " Z:" << z << " " <<
                        cc << " " << overlap << endl;
                }
            }
        }
    }
*/
    delete complete_densities;
    delete section_densities;
    delete section_reader;
    delete complete_reader;

#ifdef WITH_MPI
    MPI_Finalize();
#endif

}
