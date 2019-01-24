#include "point_transformation.h"

#include <cmath>

using std::size_t;

#define degree_cos(th) cos(M_PI/180*(th))
#define degree_sin(th) sin(M_PI/180*(th))

void point_transformation::multiply_by_matrix(double matrix[3][3], 
							double in_x, double in_y, double in_z,
							double* out_x, double* out_y, double* out_z)
{
	(*out_x) = matrix[0][0] * in_x + matrix[0][1] * in_y
									+ matrix[0][2] * in_z;
	(*out_y) = matrix[1][0] * in_x + matrix[1][1] * in_y
									+ matrix[1][2] * in_z;
	(*out_z) = matrix[2][0] * in_x + matrix[2][1] * in_y
									+ matrix[2][2] * in_z;
}

point_transformation::point_transformation(double alpha_z, double beta_y,
    double gamma_x, double translate_x,	double translate_y, double translate_z)
{
	this->set_rotation(alpha_z, beta_y, gamma_x);
	this->set_translation(translate_x, translate_y, translate_z);
}

point_transformation::point_transformation(double matrix[3][3],
    double translate_x, double translate_y, double translate_z) {
  this->set_rotation(matrix);
	this->set_translation(translate_x, translate_y, translate_z);
}

point_transformation& point_transformation::operator=(
    point_transformation other) {
  this->alpha = other.alpha;
  this->beta = other.beta;
  this->gamma = other.gamma;

  this->set_rotation(other.rotation_matrix);
  this->set_translation(other.translation[0], other.translation[1],
                        other.translation[2]);

  return *this;
}


void point_transformation::set_rotation(double alpha_z, double beta_y,
										double gamma_x)
{
	this->alpha = alpha_z;
	this->beta = beta_y;
	this->gamma = gamma_x;

	// compute cosines and sines just once
	double sin_alpha = degree_sin(alpha); 
	double sin_beta = degree_sin(beta);
	double sin_gamma = degree_sin(gamma);
	
	double cos_alpha = degree_cos(alpha); 
	double cos_beta = degree_cos(beta);
	double cos_gamma = degree_cos(gamma);
/*
	// and the negative versions for the inverse transform
	double sin_alpha_neg = degree_sin(-alpha); 
	double sin_beta_neg = degree_sin(-beta);
	double sin_gamma_neg = degree_sin(-gamma);
	
	double cos_alpha_neg = degree_cos(-alpha); 
	double cos_beta_neg = degree_cos(-beta);
	double cos_gamma_neg = degree_cos(-gamma);
*/
	// Set the 9 positions in the rotation matrix
  // Notice that since a rotation matrix is orthogonal
  // the inverse is its transpose
	rotation_matrix[0][0] = inverse_rotation_matrix[0][0]
      = cos_alpha * cos_beta;
	rotation_matrix[0][1] = inverse_rotation_matrix[1][0]
      = cos_alpha * sin_beta * sin_gamma - sin_alpha * cos_gamma;
	rotation_matrix[0][2] = inverse_rotation_matrix[2][0]
      = cos_alpha * sin_beta * cos_gamma + sin_alpha * sin_gamma;

	rotation_matrix[1][0] = inverse_rotation_matrix[0][1]
      = sin_alpha * cos_beta;
	rotation_matrix[1][1] = inverse_rotation_matrix[1][1]
      = sin_alpha * sin_beta * sin_gamma + cos_alpha * cos_gamma;
	rotation_matrix[1][2] = inverse_rotation_matrix[2][1]
      = sin_alpha * sin_beta * cos_gamma - cos_alpha * sin_gamma;
	
	rotation_matrix[2][0] = inverse_rotation_matrix[0][2]
      = - sin_beta;
	rotation_matrix[2][1] = inverse_rotation_matrix[1][2]
      = cos_beta * sin_gamma;
	rotation_matrix[2][2] = inverse_rotation_matrix[2][2]
      = cos_beta * cos_gamma;
/*
  std::cout << "Matrix:\n"
       << rotation_matrix[0][0] << "\t" << rotation_matrix[0][1] << "\t" << rotation_matrix[0][2] << "\n"
       << rotation_matrix[1][0] << "\t" << rotation_matrix[1][1] << "\t" << rotation_matrix[1][2] << "\n"
       << rotation_matrix[2][0] << "\t" << rotation_matrix[2][1] << "\t" << rotation_matrix[2][2] << "\n\n";
*/	
	// and then the inverse
/*	inverse_rotation_matrix[0][0] = cos_beta_neg * cos_alpha_neg;
	inverse_rotation_matrix[0][1] = - cos_beta_neg * sin_alpha_neg;
	inverse_rotation_matrix[0][2] = sin_beta_neg;

	inverse_rotation_matrix[1][0] = cos_gamma_neg * sin_alpha_neg
							+ sin_gamma_neg * sin_beta_neg * cos_alpha_neg;
	inverse_rotation_matrix[1][1] = cos_gamma_neg * cos_alpha_neg
							- sin_gamma_neg * sin_beta_neg * sin_alpha_neg;
	inverse_rotation_matrix[1][2] = - sin_gamma_neg * cos_beta_neg;
	
	inverse_rotation_matrix[2][0] = sin_gamma_neg * sin_alpha_neg
							- cos_gamma_neg * sin_beta_neg * cos_alpha_neg;
	inverse_rotation_matrix[2][1] = sin_gamma_neg * cos_alpha_neg
							+ cos_gamma_neg * sin_beta_neg * sin_alpha_neg;
	inverse_rotation_matrix[2][2] = cos_gamma_neg * cos_beta_neg;
*/
/*  
  std::cout << "Inverse Matrix:\n"
       << inverse_rotation_matrix[0][0] << "\t" << inverse_rotation_matrix[0][1] << "\t" << inverse_rotation_matrix[0][2] << "\n"
       << inverse_rotation_matrix[1][0] << "\t" << inverse_rotation_matrix[1][1] << "\t" << inverse_rotation_matrix[1][2] << "\n"
       << inverse_rotation_matrix[2][0] << "\t" << inverse_rotation_matrix[2][1] << "\t" << inverse_rotation_matrix[2][2] << "\n\n";
*/
}

void point_transformation::set_rotation(double matrix[3][3]) {
  for (size_t row = 0; row < 3; row++) {
    for (size_t column = 0; column < 3; column++) {
      rotation_matrix[row][column] = matrix[row][column];
      inverse_rotation_matrix[column][row] = matrix[row][column];
    }
  }
}

void point_transformation::set_translation(double translate_x,
											double translate_y,
											double translate_z)
{
	translation[0] = translate_x;
	translation[1] = translate_y;
	translation[2] = translate_z;
}

void point_transformation::transform(double in_x, double in_y,
									double in_z, double* out_x,
									double* out_y, double* out_z)
{
	double tmp_x, tmp_y, tmp_z;
	this->multiply_by_matrix(rotation_matrix, in_x, in_y, in_z,
								&tmp_x, &tmp_y, &tmp_z);
	(*out_x) = tmp_x + translation[0];
	(*out_y) = tmp_y + translation[1];
	(*out_z) = tmp_z + translation[2];
}

void point_transformation::invert_transform(double in_x, double in_y,
										double in_z, double* out_x,
										double* out_y, double* out_z)
{
	double tmp_x = in_x - translation[0];
	double tmp_y = in_y - translation[1];
	double tmp_z = in_z - translation[2];

	this->multiply_by_matrix(inverse_rotation_matrix,
								tmp_x, tmp_y, tmp_z,
								out_x, out_y, out_z);
}

string point_transformation::to_string() {
  std::stringstream ss;
  ss << "Alpha_z " << alpha << " Beta_y " << beta << " Gamma_x " << gamma
     << " tx " << translation[0] << " ty " << translation[1]
     << " tz " << translation[2];
  return ss.str();
}

void point_transformation_sequence::add(point_transformation transformation) {
  transformations.push_back(transformation); 
}

void point_transformation_sequence::add(point_transformation_sequence sequence) {
  for (size_t current = 0; current < sequence.transformations.size();
       current++) {
    this->add(sequence.transformations[current]);
  }
}
void point_transformation_sequence::add(double alpha_z, double beta_y,
    double gamma_x, double translate_x,	double translate_y,
    double translate_z) {
  point_transformation t(alpha_z, beta_y, gamma_x, translate_x, translate_y,
                         translate_z);
  this->add(t);
}
void point_transformation_sequence::add_translation(double translate_x,
    double translate_y, double translate_z) {
  this->add(0, 0, 0, translate_x, translate_y, translate_z);
}

void point_transformation_sequence::add_rotation(double alpha_z, double beta_y,
    double gamma_x) {
  this->add(alpha_z, beta_y, gamma_x, 0, 0, 0);
}

void point_transformation_sequence::clear() {
  transformations.clear();
}


void point_transformation_sequence::transform(double in_x, double in_y,
    double in_z, double* out_x, double* out_y, double* out_z) {
  if (transformations.size() == 0) {
    return;
  }
  double tmp_out_x = 0;
  double tmp_out_y = 0;
  double tmp_out_z = 0;
  double tmp_in_x = in_x;
  double tmp_in_y = in_y;
  double tmp_in_z = in_z;

  for (size_t transformation_index = 0;
       transformation_index < transformations.size(); transformation_index++) {
    point_transformation& current = transformations[transformation_index];
    current.transform(tmp_in_x, tmp_in_y, tmp_in_z,
                      &tmp_out_x, &tmp_out_y, &tmp_out_z);
    // in case there's a next iteration, update the "in" values to the
    // newly transformed ones.
    tmp_in_x = tmp_out_x;
    tmp_in_y = tmp_out_y;
    tmp_in_z = tmp_out_z;
  }

  *out_x = tmp_out_x;
  *out_y = tmp_out_y;
  *out_z = tmp_out_z;
}

void point_transformation_sequence::invert_transform(double in_x, double in_y,
    double in_z, double* out_x, double* out_y, double* out_z) {
  if (transformations.size() == 0) {
    return;
  }
  double tmp_out_x = 0;
  double tmp_out_y = 0;
  double tmp_out_z = 0;
  double tmp_in_x = in_x;
  double tmp_in_y = in_y;
  double tmp_in_z = in_z;

  for (int transformation_index = transformations.size() - 1;
       transformation_index >= 0; transformation_index--) {
    point_transformation& current = transformations[transformation_index];
    current.invert_transform(tmp_in_x, tmp_in_y, tmp_in_z,
                             &tmp_out_x, &tmp_out_y, &tmp_out_z);
    // in case there's a next iteration, update the "in" values to the
    // newly transformed ones.
    tmp_in_x = tmp_out_x;
    tmp_in_y = tmp_out_y;
    tmp_in_z = tmp_out_z;
  }

  *out_x = tmp_out_x;
  *out_y = tmp_out_y;
  *out_z = tmp_out_z;
}

point_transformation_generator::point_transformation_generator(
    double alpha_z_rotation_center, double beta_y_rotation_center,
    double gamma_x_rotation_center, double alpha_z_plus_minus,
    double beta_y_plus_minus, double gamma_x_plus_minus, double tx_center,
    double ty_center, double tz_center, double tx_plus_minus,
    double ty_plus_minus, double tz_plus_minus, double rotation_step,
    double translation_step) {
  // Initialize the first parameters that will be returned as each center
  // minus the plus_minus value for each variable.
  transformation_parameters[0] = min_values[0] = alpha_z_rotation_center -
                                                 alpha_z_plus_minus;
  transformation_parameters[1] = min_values[1] = beta_y_rotation_center -
                                                 beta_y_plus_minus;
  transformation_parameters[2] = min_values[2] = gamma_x_rotation_center -
                                                 gamma_x_plus_minus;
  transformation_parameters[3] = min_values[3] = tx_center - tx_plus_minus;
  transformation_parameters[4] = min_values[4] = ty_center - ty_plus_minus;
  transformation_parameters[5] = min_values[5] = tz_center - tz_plus_minus;
  // Max boundaries
  max_values[0] = alpha_z_rotation_center + alpha_z_plus_minus;
  max_values[1] = beta_y_rotation_center + beta_y_plus_minus;
  max_values[2] = gamma_x_rotation_center + gamma_x_plus_minus;
  max_values[3] = tx_center + tx_plus_minus;
  max_values[4] = ty_center + ty_plus_minus;
  max_values[5] = tz_center + tz_plus_minus;
  // Each variable is updated by this amount when necessary
  update_steps[0] = rotation_step;
  update_steps[1] = rotation_step;
  update_steps[2] = rotation_step;
  update_steps[3] = translation_step;
  update_steps[4] = translation_step;
  update_steps[5] = translation_step;
  // at the beginning we need to start updating the innermost loop
  next_index_to_update = 5;
}

bool point_transformation_generator::is_finished() {
  return next_index_to_update == -1;
}


point_transformation point_transformation_generator::next() {
  point_transformation current(transformation_parameters[0],
      transformation_parameters[1], transformation_parameters[2],
      transformation_parameters[3], transformation_parameters[4],
      transformation_parameters[5]);

  
  while (next_index_to_update >= 0) {
    transformation_parameters[next_index_to_update] +=
      update_steps[next_index_to_update]; 
    if (transformation_parameters[next_index_to_update] >
        max_values[next_index_to_update] ||
        min_values[next_index_to_update] ==
            max_values[next_index_to_update]) {
      // reset when we are past the max value
      // also, if the min and max are equal we just ignore it and continue
      transformation_parameters[next_index_to_update] =
          min_values[next_index_to_update];
      // go to the outer loop
      next_index_to_update--;
    } else { // the update didn't end the current loop
      // reset the inner loops to the min value and start from the innermost
      while (next_index_to_update < 5) {
        next_index_to_update++;
        transformation_parameters[next_index_to_update] =
            min_values[next_index_to_update];
      }
      break;
    }
  }
  return current;
}
/*
#include <iostream>

using std::cout;

int main() {
  point_transformation_generator g(180,0,180, // centers
      180, 90, 180, 0,0,0, 0, 0, 0, 90, 0);
  while(!g.is_finished()) {
    point_transformation t = g.next();
    cout << t.to_string() << "\n";
  }
*/
/*
  for (double alpha_z = -30; alpha_z < 30; alpha_z+= 5) {
    for (double beta_y = -15; beta_y < 15; beta_y+= 5) {
      for (double gamma_x = -10; gamma_x < 10; gamma_x += 5) {
        for (double x = -5; x < 5; x++) {
          for (double y = -3; y < 3; y++) {
            for (double z = -2; z < 2; z++) {
              point_transformation t(alpha_z, beta_y, gamma_x, x, y, z);
              cout << t.to_string() << "\n";
            }
          }
        }
      }
    }
  }
*/
/*  return 0;
}*/

/*
#include <iostream>

using std::cout;

int main()
{
  point_transformation t1(5.22,15.81,25.22,12.4533,87.43,-12.554);
  point_transformation t2(10,63,55,-87,-87.2456,-100.3442);

  point_transformation_sequence seq;
  seq.add(t1);
  seq.add(t2);

  double p1[3] = {55.37, -5.91, 12.46};
	double out[3];
  double inv[3];
	std::cout << p1[0] << " " << p1[1] << " " << p1[2] << "\n";
  
	seq.transform(p1[0], p1[1], p1[2], &out[0], &out[1], &out[2]);
	std::cout << out[0] << " " << out[1] << " " << out[2] << "\n";
	
  seq.invert_transform(out[0], out[1], out[2], &inv[0], &inv[1], &inv[2]);
	std::cout << inv[0] << " " << inv[1] << " " << inv[2] << "\n";

  cout << "-----------------\n";

	t1.transform(p1[0], p1[1], p1[2], &out[0], &out[1], &out[2]);
	std::cout << out[0] << " " << out[1] << " " << out[2] << "\n";
  
  t2.transform(out[0], out[1], out[2], &inv[0], &inv[1], &inv[2]);
	std::cout << inv[0] << " " << inv[1] << " " << inv[2] << "\n";
}


int main()
{
//	point_transformation t(-53.25, -114.38, -138.05, 0, 0, 0);
//	point_transformation t(41.94505915,65.61575591, 233.2517417, 0, 0, 0);
//  point_transformation t(156.0078714, -146.3543482, 10.41251614,0,0,0);
//  point_transformation t(169.5874839, 33.64565176, 23.9921286, 0, 0, 0);
//  point_transformation t(23.9921286,33.64565176, 169.5874839, 0, 0, 0);
  point_transformation t(15,15,15,0,0,0);

//	double p1[3] = {18.625, 51.764, -14.042};
  double p1[3] = {1.0f, 1.0f, 1.0f};
	double out[3];
//	double inv[3];

	t.transform(p1[0], p1[1], p1[2], &out[0], &out[1], &out[2]);
	std::cout << out[0] << " " << out[1] << " " << out[2] << "\n";
  p1[0]++;p1[1]++;p1[2]++;

	t.transform(p1[0], p1[1], p1[2], &out[0], &out[1], &out[2]);
	std::cout << out[0] << " " << out[1] << " " << out[2] << "\n";
  p1[0]++;p1[1]++;p1[2]++;

	t.transform(p1[0], p1[1], p1[2], &out[0], &out[1], &out[2]);
	std::cout << out[0] << " " << out[1] << " " << out[2] << "\n";
  p1[0]++;p1[1]++;p1[2]++;

	t.transform(p1[0], p1[1], p1[2], &out[0], &out[1], &out[2]);
	std::cout << out[0] << " " << out[1] << " " << out[2] << "\n";
  p1[0]++;p1[1]++;p1[2]++;

	t.transform(p1[0], p1[1], p1[2], &out[0], &out[1], &out[2]);
	std::cout << out[0] << " " << out[1] << " " << out[2] << "\n";
  p1[0]++;p1[1]++;p1[2]++;

	t.transform(p1[0], p1[1], p1[2], &out[0], &out[1], &out[2]);
	std::cout << out[0] << " " << out[1] << " " << out[2] << "\n";
  p1[0]++;p1[1]++;p1[2]++;
//	t.invert_transform(out[0], out[1], out[2], &inv[0], &inv[1], &inv[2]);
	
//	std::cout << inv[0] << " " << inv[1] << " " << inv[2] << "\n";

}
*/
