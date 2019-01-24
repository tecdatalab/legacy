#include "density_lattice.h"
#include "../proteinstructure/point_transformation.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>

#ifdef WITH_MPI

// only include if the mpi version is being compiled
#include "mpi.h"

#endif

using namespace std;

double density_lattice::calculate_correlation(density_lattice& compared_lattice,
    point_transformation_sequence& compared_transformation,
    density_vector min_bound, density_vector max_bound, float density_threshold,
    double* overlap) {
  float this_density, compared_density;
	unsigned long voxel_overlap = 0;
  unsigned long compared_lattice_voxels_over_threshold = 0;
  double sum_of_products = 0;
  double this_sum_squared = 0;
  double compared_sum_squared = 0;
  double tx, ty, tz;

  for (double x = min_bound.x; x <= max_bound.x; x += voxel_length_x()) {
    for (double y = min_bound.y; y <= max_bound.y; y += voxel_length_y()) {
      for (double z = min_bound.z; z <= max_bound.z; z += voxel_length_z()) {
        compared_transformation.invert_transform(x, y, z, &tx, &ty, &tz);
        this_density = this->get_from_real_coordinate(x, y, z);
        compared_density =
            compared_lattice.get_from_real_coordinate(tx, ty, tz);
        if(compared_density >= density_threshold) {
          compared_lattice_voxels_over_threshold++;
				  if (this_density >= density_threshold) {
            sum_of_products += this_density * compared_density;
            this_sum_squared += this_density * this_density;
            compared_sum_squared += compared_density * compared_density;
            voxel_overlap++;
          }
        }
      }
    }
  }
	// No overlap => assign neutral values as result
	if (voxel_overlap == 0) {
		*overlap = 0;
		return 0;
	} else {
		// set the overlap fraction
		*overlap = static_cast<double>(voxel_overlap) /
               compared_lattice_voxels_over_threshold;
    double denominator = sqrt(this_sum_squared * compared_sum_squared);
    return sum_of_products / denominator;
  }
}

// TODO: This second version is legacy code. When sampling is clearly defined
// merge the two or remove the second one.
double density_lattice::calculate_correlation(density_lattice* compared_lattice,
                                              unsigned long offset_x,
                                              unsigned long offset_y,
                                              unsigned long offset_z,
								 			                        double alpha_z, double beta_y,
			                      								  double gamma_x,
                                              float density_threshold,
                                              unsigned long* overlap)
{
  float this_density, compared_density;
	unsigned long voxel_overlap = 0;
  double sum_of_products = 0;
  double this_sum_squared = 0;
  double compared_sum_squared = 0;
  double tx, ty, tz;

	// this is used to determine the compared_lattice point required,
	// given that we are translating and rotating it
	point_transformation transformation(alpha_z, beta_y, gamma_x,
                  										offset_x, offset_y, offset_z);
    
  for (unsigned long x = 0; x < this->size_x(); x++) {
    for (unsigned long y = 0; y < this->size_y(); y++) {
      for (unsigned long z = 0; z < this->size_z(); z++) {
				transformation.transform(x, y, z, &tx, &ty, &tz);
        compared_density = compared_lattice->get_interpolated(tx, ty, tz);
				// and get the density in this lattice for comparison
        this_density = this->get(x, y, z);
				if (this_density >= density_threshold &&
            compared_density >= density_threshold) {
          sum_of_products += this_density * compared_density;
          this_sum_squared += this_density * this_density;
          compared_sum_squared += compared_density * compared_density;
          voxel_overlap++;
        }
      }
    }
  }
	// No overlap => assign neutral values as result
	if (voxel_overlap == 0) {
		*overlap = 0;
		return 0;
	} else {
		// set the overlap fraction
		*overlap = voxel_overlap;
    double denominator = sqrt(this_sum_squared * compared_sum_squared);
    return sum_of_products / denominator;
  }
}

void density_lattice::get_min_max_densities(float* min,
    float* max) {
  unsigned long x, y, z;
  float partial_min = numeric_limits<double>::max();
  float partial_max = -numeric_limits<double>::max();
  this->iterator_reset();
  while (!this->iterator_is_done()) {
    float density = this->iterator_get_next_value(&x, &y, &z);
    if (density < partial_min) {
      partial_min = density;
    }
    if (density > partial_max) {
      partial_max = density;
    }
  }
  *min = partial_min;
  *max = partial_max;
}

void density_lattice::get_voxel_boundaries_for_threshold(
    float density_threshold, density_voxel* min_bound,
    density_voxel* max_bound) {
	float density;
  unsigned long x, y, z, min_x, min_y, min_z, max_x, max_y, max_z;
	// Assign the boundaries to edge values
	min_x = min_y = min_z = ULONG_MAX;
	max_x = max_y = max_z = 0;
	
  this->iterator_reset();
  while (!this->iterator_is_done()) {
    density = this->iterator_get_next_value(&x, &y, &z);
    if (density > density_threshold) {
      if (x < min_x) min_x = x;
      if (y < min_y) min_y = y;
      if (z < min_z) min_z = z;
            
			if (x > max_x) max_x = x;
      if (y > max_y) max_y = y;
      if (z > max_z) max_z = z;
    }
  }

  (*min_bound).x = min_x;
  (*min_bound).y = min_y;
  (*min_bound).z = min_z;
	(*min_bound).density = this->get(min_x, min_y, min_z);
  
  (*max_bound).x = max_x;
  (*max_bound).y = max_y;
  (*max_bound).z = max_z;
	(*max_bound).density = this->get(max_x, max_y, max_z);
}

void density_lattice::get_boundaries_for_threshold(float density_threshold,
    density_vector* min_bound, density_vector* max_bound) {
  density_voxel min, max;
  get_voxel_boundaries_for_threshold(density_threshold, &min, &max);

//cout << "VOXELS " << min.x << " " << min.y << " " << min.z << "\n";
  // Adjust the voxel numbers to actual coordinates and set the densities.
  (*min_bound).x = (static_cast<double>(min.x) + origin_x()) * voxel_length_x();
  (*min_bound).y = (static_cast<double>(min.y) + origin_y()) * voxel_length_y();
  (*min_bound).z = (static_cast<double>(min.z) + origin_z()) * voxel_length_z();
	(*min_bound).density = min.density;

  (*max_bound).x = (static_cast<double>(max.x) + origin_x()) * voxel_length_x();
  (*max_bound).y = (static_cast<double>(max.y) + origin_y()) * voxel_length_y();
  (*max_bound).z = (static_cast<double>(max.z) + origin_z()) * voxel_length_z();
	(*max_bound).density = max.density;
}

vector<density_vector> density_lattice::calculate_all_corners(
    const density_vector& min, const density_vector& max) {
  vector<density_vector> corners;
  // Temporary holder of coordinates. 6 copies of it will be pushed
  // into the vector after the various modifications to get each corner.
  density_vector current;
  // The parameters are the two trivial corners we already have
  corners.push_back(min);
  corners.push_back(max);
  // Start from the min coordinates and change each value to create 3 new
  // corners, one per axis displacement
  current = min;
  current.x = max.x;
  corners.push_back(current);
  current = min;
  current.y = max.y;
  corners.push_back(current);
  current = min;
  current.z = max.z;
  corners.push_back(current);
  // Do the same starting from the max, replacing each axis value by min.?
  current = max;
  current.x = min.x;
  corners.push_back(current);
  current = max;
  current.y = min.y;
  corners.push_back(current);
  current = max;
  current.z = min.z;
  corners.push_back(current);

  return corners;
}


void density_lattice::get_boundaries_after_transformation(
    point_transformation_sequence transform,
    density_vector min_before, density_vector max_before,
    density_vector* min_after, density_vector* max_after) {
  // 1. Get all corners before transformation.
  vector<density_vector> corners = calculate_all_corners(min_before,
                                                         max_before);
  // 2. Transform them and scan each, keeping the lowest x, y and z
  // values.
  double transformed_x, transformed_y, transformed_z, min_x, min_y, min_z,
         max_x, max_y, max_z;
  min_x = min_y = min_z = numeric_limits<double>::max();
  max_x = max_y = max_z = -numeric_limits<double>::max();


  for (size_t corner_index = 0; corner_index < corners.size(); corner_index++) {
    const density_vector& current_corner = corners[corner_index];
//cout << "Corner " <<current_corner.x << " " << current_corner.y
//     << " " << current_corner.z << "\n";

    transform.transform(current_corner.x, current_corner.y, current_corner.z,
                        &transformed_x, &transformed_y, &transformed_z);
//cout << "T " <<transformed_x << " " << transformed_y
//     << " " << transformed_z << "\n";
    if (transformed_x < min_x) min_x = transformed_x;
    if (transformed_y < min_y) min_y = transformed_y;
    if (transformed_z < min_z) min_z = transformed_z;
    if (transformed_x > max_x) max_x = transformed_x;
    if (transformed_y > max_y) max_y = transformed_y;
    if (transformed_z > max_z) max_z = transformed_z;
  }
  // Assign the final values to the output parameters.
  min_after->x = min_x;
  min_after->y = min_y;
  min_after->z = min_z;
  max_after->x = max_x;
  max_after->y = max_y;
  max_after->z = max_z;
/*
cout << "Final " << min_x << " " << min_y << " " << min_z
     << " " << max_x << " " << max_y << " " << max_z << "\n";*/
}

correlation_results density_lattice::get_significant_correlations(
	density_lattice* compared_lattice, int step_size, double rotation_step,
	float density_threshold, float correlation_threshold,
	float overlap_threshold)
{
	correlation_results all_results;
	correlation_result current_result;
	// The translations start and finish at the volume boundaries for
	// the threshold provided
	density_voxel min_bound, max_bound, this_min_bound, this_max_bound,
                compared_min_bound, compared_max_bound;
  // This represents the denominator in the overlap fraction calculation
  float voxels_over_threshold = count_voxels_over_threshold(density_threshold);
	compared_lattice->get_voxel_boundaries_for_threshold(density_threshold,
                                                       &compared_min_bound,
                                                       &compared_max_bound);
	this->get_voxel_boundaries_for_threshold(density_threshold, &this_min_bound,
                                           &this_max_bound);
  min_bound.x = max(0ul, compared_min_bound.x - this_min_bound.x);
  min_bound.y = max(0ul, compared_min_bound.y - this_min_bound.y);
  min_bound.z = max(0ul, compared_min_bound.z - this_min_bound.z);
  max_bound = compared_max_bound; 
  // Calculate the total number of CC computations that would be done.
  // For non-MPI runs this is just used to output progress.
  // For MPI Have all the ranks go through all x,y,z,alpha,beta,gamma values
  // but have them ignore all except the ones that are in their range
  long long total_computations = 
    static_cast<long long>(
      ceil((max_bound.x - min_bound.x + 1) / float(step_size)) *
      ceil((max_bound.y - min_bound.y + 1) / float(step_size)) *
      ceil((max_bound.z - min_bound.z + 1) / float(step_size)) *
      ceil(360 / rotation_step) *
      ceil(360 / rotation_step) *
      ceil(180 / rotation_step));

  if (total_computations < 0) {
    cerr << "Overflow in comparisons counter, density_lattice.get_significant_correlations" << endl;
  }

#ifdef WITH_MPI
  // If mpi is used the different bases will be distributed between
  // different ranks
  int int_rank, int_numtasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &int_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &int_numtasks);
  // Have long long variables to force 64 bits. The number of total_computations
  // can be really large and these two variables are used in related calculations.
  long long rank = static_cast<long long>(int_rank);
  long long numtasks = static_cast<long long>(int_numtasks);
  long long iterations_per_task = total_computations / numtasks;
  long long total_remainder = total_computations % numtasks;
  long long remainder_processed = total_remainder > rank ? rank : total_remainder;
  long long start_index = iterations_per_task * rank + remainder_processed;
  long long end_index = start_index + iterations_per_task;
  // the remainder will be given 1 by 1 to each of the first
  // ranks, until we run out of remaining bases
  if(rank < total_remainder)
  {
	  end_index++;
  }

#endif
  // Keep track of the iteration number with this variable
  long long current_iteration = 0;
  // This one is used to print out on 1% intervals
  double completion_percentage = 0.0;
  // This one is calculated every iteration
  double exact_completion_percentage = 0.0;

	// 6D exhaustive search
	for (unsigned long x = min_bound.x; x <= max_bound.x; x+=step_size)
	{
  	for (unsigned long y = min_bound.y; y <= max_bound.y; y+=step_size)
		{
    	for (unsigned long z = min_bound.z; z <= max_bound.z; z+=step_size)
			{
				for (double alpha_z = 0; alpha_z < 360; alpha_z+=rotation_step)
				{
					for (double beta_y = -90; beta_y < 90; beta_y+=rotation_step)
					{
						for (double gamma_x = 0; gamma_x < 360;
								gamma_x += rotation_step)
						{
              current_iteration++;
#ifdef WITH_MPI
              if ((current_iteration - 1) < start_index ||
                  (current_iteration - 1) >= end_index) {
                // Ignore it, this rank is not supposed to compute it
                continue;
              }
              exact_completion_percentage = (current_iteration - start_index) /
                  static_cast<double>(end_index - start_index);
              // report every 1%
              if (exact_completion_percentage - completion_percentage > 0.01) {
                completion_percentage += 0.01;
                cerr << "Process " << rank << " " << (completion_percentage * 100)
                     << "\% complete\n";
              }
#else
              exact_completion_percentage = current_iteration /
                  static_cast<double>(total_computations);
              // report every 1%
              if (exact_completion_percentage - completion_percentage > 0.01) {
                completion_percentage += 0.01;
                cerr << (completion_percentage * 100) << "\% complete\n";
              }
#endif
					    unsigned long overlap = 0;
							double cc = this->calculate_correlation(compared_lattice,	x, y, z,
                                                      alpha_z, beta_y, gamma_x,
																				              density_threshold,
          																						&overlap);
              float overlap_fraction = float(overlap) / voxels_over_threshold;
							if (cc >= correlation_threshold &&
				          overlap_fraction >= overlap_threshold) {
								current_result.x = x;
								current_result.y = y;
								current_result.z = z;
								current_result.alpha_z = alpha_z;
								current_result.beta_y = beta_y;
								current_result.gamma_x = gamma_x;
								current_result.correlation = cc;
								current_result.overlap = overlap_fraction;
								all_results.push_back(current_result);
              }
					  }
					}
				}
      }
    }
  }
#ifdef WITH_MPI
  cerr << "Process " << rank << " 100\% complete\n";
#else
  cerr << "100\% complete\n";
#endif
	return all_results;
}

double density_lattice::calculate_volume(float density_threshold)
{
	double voxel_volume = this->voxel_length_x() *
						  this->voxel_length_y() *
						  this->voxel_length_z();
  double voxel_count = this->count_voxels_over_threshold(density_threshold);
  return voxel_volume * voxel_count;
}

unsigned long density_lattice::count_voxels_over_threshold(
										float density_threshold)
{
#ifdef DEBUG
cerr << "Counting voxels threshold " << density_threshold << endl;
#endif
	unsigned long count = 0;
	iterator_reset();
	while (!iterator_is_done()) {
		if(iterator_get_next_value() >= density_threshold) {
			count++;
		}
	}
	return count;
}
