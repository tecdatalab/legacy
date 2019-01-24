#include "main_em_correlation_clustering_options.h"

#include "density_lattice.h"
// Both .h and .cc are kept for k-means for standardization purposes, however,
// since they contain templates and they are expressed in this file, both
// .h and .cc must be included/compiled as part of this main program.
#include "kmeans.h"
#include "kmeans.cc"
#include "mrc.h"

#include <cstdlib>
#include <iostream>

/*
 * Clusters results from em_correlation by computing the x,y,z coordinates
 * of the EM centers (after transformation) and clustering them using k-means
 */
int main(int argc, char** argv)
{
	main_em_correlation_clustering_options options(argc,argv);
  if(!options.parse_successful()) {
  	cout << options.usage();
    return 0;
  }
  // 1. Calculate the center of the EM map volume.
  // These are the x,y,z coordinates that will be transformed to
  // generated the potential cluster centers.
  mrc_reader* em_reader = new mrc_reader(options.get_em_filename());
  density_lattice* em_lattice = em_reader->get_densities();
  float density_threshold = options.get_contour_value();
  density_vector min_bound, max_bound;

  em_lattice->get_boundaries_for_threshold(density_threshold,
                                           &min_bound, &max_bound);
  double center_x = ((max_bound.x - min_bound.x) / 2.0f) + min_bound.x;
  double center_y = ((max_bound.y - min_bound.y) / 2.0f) + min_bound.y;
  double center_z = ((max_bound.z - min_bound.z) / 2.0f) + min_bound.z;

  delete em_lattice;
  delete em_reader;

  // 2. Load the existing results, calculate where the center would lie
  // given each transformation, and cluster them.
  correlation_results input_results;
  input_results.read(options.get_input_filename());
  kmeans<correlation_result> cluster_manager(3); // 3 dimensions (x,y,z)

  for(size_t result_index = 0; result_index < input_results.size();
      result_index++) {
    input_results[result_index].calculate_transformed_center(center_x,
                                                             center_y,
                                                             center_z);
    cluster_manager.add_node(&input_results[result_index]);
  }
  // 2.1. Calculate the silhouette coefficient for each value of k.
  size_t best_k = 1;
  double best_coefficient = -1; // -1 is the worst possible value.

  cout << "k\tSilhouette\tRadius\tDiameter" << endl;

  for(size_t k = 2; k <= options.get_k(); k++) {
    cluster_manager.initialize(k);
    cluster_manager.cluster();

    double current_coefficient =
      cluster_manager.calculate_silhouette_coefficient();

    if(current_coefficient > best_coefficient) {
      best_k = k;
      best_coefficient = current_coefficient;
    }
    cout << k << "\t" << current_coefficient 
         << "\t" << cluster_manager.calculate_average_cluster_radius()
         << "\t" << cluster_manager.calculate_average_cluster_diameter()
         << endl;
  }

  // 2.2. Re-cluster using the best k found
  cluster_manager.initialize(best_k);
  cluster_manager.cluster();

  // 3. Output the cluster centers.
  vector<correlation_result*> centroids = cluster_manager.get_centroids();
  correlation_results clustered_results;
  for(size_t centroid_index = 0; centroid_index < centroids.size();
      centroid_index++) {
    clustered_results.push_back(*centroids[centroid_index]);
  }

  clustered_results.write(options.get_output_filename());
}
