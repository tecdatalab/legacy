#include "density_lattice.h"
#include "../proteinstructure/point_transformation.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <list>

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
	return this->count_voxels_over_threshold(density_threshold) *
			voxel_volume;
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



/*
 * Start of Additional Code to Remove Noise 
 * from EM density map
 *
 * Class Merge_Manager
 * Responsible for maintaining a list of sets.
 * Each set holds group numbers of those groups that touch
 * in the density lattice object.
 *
 */

class merge_manager {
	list< std::set< unsigned int > > merge_groups;

	//Returns true if an element in 'add' is also in 'present'
	bool intersects(const set<unsigned int>& add, const set<unsigned int>& present) {
		for (set<unsigned int>::iterator itr=add.begin(); itr != add.end(); ++itr) {
			if (present.find(*itr) != present.end())
				return true;
		}
		return false;
	}

public:
	/*
	 * Check each set in the list merge_groups. If an element (group number)
	 * in add_set matches an element in the list's set, remove the set
	 * from merge_groups and add it into the add_set. Do this for
	 * all sets in the list merge_groups. Once completed, add
	 * add_set into the list merge_groups.
	 */
	void add_set(std::set<unsigned int>& add_set) {
		//itr on list, points to set(s)
		list< std::set< unsigned int > >::iterator itr = merge_groups.begin();
		while (itr != merge_groups.end()){
			if (intersects(add_set, (*itr))) {
				//insert needs 2 set itr's, pointing to uint(s)
				add_set.insert(itr->begin(), itr->end());
				itr = merge_groups.erase(itr);
			}
			else {//no shared group number found; check next set in the list
				++itr;
			}
		}
		//add add_set into merge_groups
		merge_groups.push_front(add_set);
	}

	/*
	 * Returns a pointer to the set that contains the given int;
	 * otherwise returns NULL.
	 */
	std::set<unsigned int>* contains(unsigned int index) {

		for (list< std::set< unsigned int > >::iterator itr = merge_groups.begin(); itr != merge_groups.end(); itr++) {
			if (itr->find(index) != itr->end()) {
				return &(*itr);
			}
		}
		return NULL;
	}
};

/*
 * Main function for removing noise from EM density map.
 * Creates groups matrix; calls "assign_group" for each
 * (x, y, z) coordinate; calls delete_groups on the 
 * minimum group numbers. 
 */
void density_lattice::remove_noise(float minimum) {


	//count[index #] == # of (x, y, z) coordinates that are in group/index #
	std::map<unsigned int, unsigned int> counts;

	merge_manager merge_groups;

	//Current group to assign next (x,y,z) value into
	unsigned int next_group = 1;

	long unsigned int x_size =  this->size_x();
	long unsigned int y_size = this->size_y();
	long unsigned int z_size = this->size_z();
	//Create groups matrix: same size as values, initialized to zero
	unsigned int ***groups;

	groups = new unsigned int**[x_size];
	for (unsigned long x = 0; x < x_size; x++) {
		groups[x] = new unsigned int*[y_size];
		for (unsigned long y = 0; y < y_size; y++) {
			// the parentheses trigger the default initialization (zeroes)
			groups[x][y] = new unsigned int[z_size]();
		}
	}

	//Check 6 neighbors of each (x, y, z) to assign groups
	for (unsigned long x = 0; x < x_size; x++) {
		for (unsigned long y = 0; y < y_size; y++) {
			for (unsigned long z = 0; z < z_size; z++) {
				assign_group(x, y, z, groups, next_group, counts, minimum, merge_groups);
			}
		}
	}

	cout << "Finished assigning groups; beginning deletion\n";

	//Keep dominant groups by assigning smaller groups to have min value
	delete_small_groups(groups, minimum, counts, merge_groups);

	//Free up groups
	for (unsigned long x = 0; x < x_size; x++) {
		for (unsigned long y = 0; y < y_size; y++) {
			delete [] groups[x][y];
		}
		delete [] groups[x];
	}
	delete [] groups; 
}

/*
 * Calculates the number of values within each group.
 * Saves the dominant group by overwriting members of other 
 * groups to have the dmin value in the values matrix.
 */
void density_lattice::delete_small_groups(
		unsigned int ***groups, float min,
		map<unsigned int, unsigned int>& count, merge_manager& merge_groups) {

	//Something went wrong 
	if (count.size() < 2) {
		cout << "Error in count hashmap\n";
		return;
	}


	unsigned int max_value = 0;
	unsigned int max_group = 0;
	std::set<unsigned int> *max_groups;

	//Check the map of counts starting with group #1 (skipping 0 = min group)
	for (map<unsigned int, unsigned int>::const_iterator itr = ++count.begin(); itr != count.end(); itr++) {
		if (itr->second != 0) {
			std::set<unsigned int> *group = merge_groups.contains(itr->first);
			if (group != NULL) {
				std::set<unsigned int>::iterator first_itr = group->begin();
				for (std::set<unsigned int>::iterator set_itr = ++group->begin(); set_itr != group->end(); set_itr++) {
					count[*first_itr] += count[*set_itr];
					count[*set_itr] = 0;
				}
				if (count[*first_itr] > max_value) {
					max_value = count[*first_itr];
					max_groups = group;
					max_group = 0;
				}
			}
			//if group was null, then it shares no group members
			if (itr->second > max_value) {
				max_value = itr->second;
				max_group = itr->first;
				max_groups = NULL;
			}
		}
	}
	
	//Case when there's only one large, isolated group
cout << "Max group ID " << max_group << endl;
	if (max_group != 0) {
		cout << "Setting groups using max_group: " << max_group << "\n";
		//Set all values in the other groups to be the minimum value (min)
		for (unsigned long x = 0; x < this->size_x(); x++) {
			for (unsigned long y = 0; y < this->size_y(); y++) {
				for (unsigned long z = 0; z < this->size_z(); z++) {
					if (groups[x][y][z] != max_group)
						this->set(min, x, y, z);
					//values[x][y][z] = min;
				}
			}
		}
	}
	
	//Case when there are multiple groups touching
	if (max_groups != NULL) {
		cout << "Setting groups using max_groups\n";
		//Set all values in the other groups to be the minimum value (min)
		for (unsigned long x = 0; x < this->size_x(); x++) {
			for (unsigned long y = 0; y < this->size_y(); y++) {
				for (unsigned long z = 0; z < this->size_z(); z++) {
					if (max_groups->find(groups[x][y][z]) == max_groups->end())
						this->set(min, x, y, z);
				}
			}
		}
	}
}


/*
 * Assigns a group number to the given (x, y, z) coordinate
 * based on the group numbers of its 6 surrounding neighbors.
 */
void density_lattice::assign_group(unsigned long x,
		unsigned long y, unsigned long z, unsigned int ***groups,
		unsigned int& next_group, std::map<unsigned int, unsigned int>& count,
		float minimum, merge_manager& merge_groups) {

	//cout << "Value at x: " << x << ", y: " << y << ", z: " << z << " is: " << this->get(x,y,z) << "\n";
	//cout << "Dmin is: " << minimum << "\n";

	//2 general cases: min vs. not min
	if (this->get(x,y,z) != minimum){

		//6 neighbors: 0 - backward, 1 - forward, 2 - down, 3 - up, 4 - left, 5 - right
		std::set<unsigned int> neighbor_groups = check_neighbors(x, y, z, groups);

		if (neighbor_groups.empty()) { //all neighbors are group 0 (min value)
			//cout << "No neighbors; using next_group: " << next_group << "\n";
			groups[x][y][z] = next_group;
			count[next_group] = 1;
			next_group++;
		}
		else { //at least one neighbor is in a non-zero group

			unsigned int group = *(neighbor_groups.begin());
			//cout << "Some neighbors; using group: " << group << "\n";

			groups[x][y][z] = group;
			count[group]++;

			if (neighbor_groups.size() > 1)
				merge_groups.add_set(neighbor_groups);
		}
	}
	else { //(values[x][y][z] == dmin -> minimum values into the zero group
		//cout << "Assigning x: " << x << ", y: " << ", z: " << z << " to group 0 \n";
		groups[x][y][z] = 0;
		count[0]++; //not necessary to track this really
	}
}


/*
 * Checks the neighboring coordinates to determine what the new
 * current group number should be.
 * For instance if there is a pocket of non-min values surrounded
 * by min values (the zero group), but there has already been a
 * series of non-min values placed into another group, say group 1,
 * then this pocket should be assigned to group 2, under the condition
 * that no members of this pocket touch any members of group 1.
 * If they do touch, then members of this pocket should be assigned to
 * the next_group.
 *
 * Know there must be non-zero group neighbors.
 *

list<unsigned int> density_lattice::get_neighboring_groups(
		unsigned long x, unsigned long y, unsigned long z,
		unsigned int* neighbor_groups) {

	int i,j = 0;
	unsigned int index = 0;
	unsigned int seen[6];

	for (i = 0; i < 6; i++) {
		seen[i] = 0;

		if (neighbor_groups[i] == 0) continue;

		for (j = i - 1; j >= 0; j--) {
			if (neighbor_groups[i] == neighbor_groups[j]) {
				seen[i] = seen[j];
				break;
			}
		}
		seen[i]++;
		if (seen[i] > seen[index])
			index = i;
	}

	cout << "Returning new group number: " << neighbor_groups[index] << "\n";
	return neighbor_groups[index];
}
 */

/*
 * Fills the "neighbor_groups" set with the group number
 * of the 6 neighbors of coordinate (x,y,z).
 *
 * Previous implementation: A value of zero in "neighbor_groups" indicates
 * no neighbor in that direction, i.e. a boundary case.
 */
set<unsigned int> density_lattice::check_neighbors(unsigned long x,
		unsigned long y, unsigned long z,
		unsigned int ***groups) {

	std::set<unsigned int> neighbor_groups;

	if(x != 0 && groups[x-1][y][z] != 0 )
		neighbor_groups.insert(groups[x-1][y][z]);
	if (x != (this->size_x()-1) && groups[x+1][y][z] != 0)
		neighbor_groups.insert(groups[x+1][y][z]);
	if (y != 0 && groups[x][y-1][z] != 0)
		neighbor_groups.insert(groups[x][y-1][z]);
	if (y != (this->size_y()-1) && groups[x][y+1][z] != 0 )
		neighbor_groups.insert(groups[x][y+1][z]);
	if (z != 0 && groups[x][y][z-1] != 0)
		neighbor_groups.insert(groups[x][y][z-1]);
	if (z != (this->size_z()-1) &&  groups[x][y][z+1] != 0)
		neighbor_groups.insert(groups[x][y][z+1]);

	return neighbor_groups;
}
