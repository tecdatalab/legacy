#include "point_lattice.h"

#include <cmath>
#include <iostream>

point_lattice::point_lattice(unsigned long in_nx, unsigned long in_ny,
    unsigned long in_nz, float in_xorigin, float in_yorigin, float in_zorigin,
    float in_voxel_xlen, float in_voxel_ylen, float in_voxel_zlen) :
        nx(in_nx), ny(in_ny), nz(in_nz),
        xorigin(in_xorigin), yorigin(in_yorigin), zorigin(in_zorigin),
        voxel_xlen(in_voxel_xlen), voxel_ylen(in_voxel_ylen),
        voxel_zlen(in_voxel_zlen), iteratorx(0), iteratory(0), iteratorz(0)
{
    values = new float**[nx];
    for (unsigned long x = 0; x < nx; x++) {
        values[x] = new float*[ny];
        for (unsigned long y = 0; y < ny; y++) {
            // the parentheses trigger the default initialization (zeroes)
            values[x][y] = new float[nz]();
        }
    }
}

point_lattice::~point_lattice()
{
    for (unsigned long x = 0; x < nx; x++) {
        for (unsigned long y = 0; y < ny; y++) {
            delete[] values[x][y];
        }
        delete[] values[x];
    }
    delete[] values;
}


void point_lattice::real_coordinate_to_voxel(double real_x, double real_y,
    double real_z, unsigned long* voxel_x, unsigned long* voxel_y,
    unsigned long* voxel_z) {
  // +0.5f to round to the closest integer
  signed long raw_x = static_cast<unsigned long>(
      floor((real_x - (xorigin * voxel_xlen)) / voxel_xlen + 0.5f));
  signed long raw_y  = static_cast<unsigned long>(
      floor((real_y - (yorigin * voxel_ylen)) / voxel_ylen + 0.5f));
  signed long raw_z = static_cast<unsigned long>(
        floor((real_z - (zorigin * voxel_zlen)) / voxel_zlen + 0.5f));
  *voxel_x = raw_x >= 0 && static_cast<unsigned long>(raw_x) < nx ?
                static_cast<unsigned long>(raw_x) :
                raw_x < 0 ? 0 : nx -1;
  *voxel_y = raw_y >= 0 && static_cast<unsigned long>(raw_y) < ny ?
                static_cast<unsigned long>(raw_y) :
                raw_y < 0 ? 0 : ny -1;
  *voxel_z = raw_z >= 0 && static_cast<unsigned long>(raw_z) < nz ?
                static_cast<unsigned long>(raw_z) :
                raw_z < 0 ? 0 : nz -1;
}

void point_lattice::voxel_to_real_coordinate(unsigned long voxel_x,
    unsigned long voxel_y, unsigned long voxel_z, double* real_x,
    double* real_y, double* real_z) {
  *real_x = voxel_x * voxel_xlen + xorigin * voxel_xlen;
  *real_y = voxel_y * voxel_ylen + yorigin * voxel_ylen;
  *real_z = voxel_z * voxel_zlen + zorigin * voxel_zlen;
}

// Internally calls get_interpolated but it first scales down the real
// coordinate and translates it to make it a coordinate in voxel space.
float point_lattice::get_from_real_coordinate(double x, double y, double z) {
  float voxel_x = (x - (xorigin * voxel_xlen)) / voxel_xlen;
  float voxel_y = (y - (yorigin * voxel_ylen)) / voxel_ylen;
  float voxel_z = (z - (zorigin * voxel_zlen)) / voxel_zlen;
  return get_interpolated(voxel_x, voxel_y, voxel_z);
}

float point_lattice::get_interpolated(float x, float y, float z)
{
	// Note: the notation here follows the trilinear interpolation
	// explanation found on wikipedia, where x0, y0 and z0 represent
	// the lowest x, y, z coordinates in the voxel and xd, yd, zd
	// the offset that needs to be added along each axis
	// the get to the interpolated coordinates
	unsigned long x0 = (unsigned long)floor(x);
	unsigned long y0 = (unsigned long)floor(y);
	unsigned long z0 = (unsigned long)floor(z);
	//std::cerr << x << " " << y << " " << z << std::endl;
	// If out of bounds we'll just return zero
	if(x0 >= (nx - 1) || y0 >= (ny - 1) || z0 >= (nz - 1) ||
		floor(x) < 0 || floor(y) < 0 || floor(z) < 0) {
		return 0.0;
	}

	unsigned long x1 = x0 + 1;
	unsigned long y1 = y0 + 1;
	unsigned long z1 = z0 + 1;

	float xd = (x - (float)x0) / (x1 - x0);
	float yd = (y - (float)y0) / (y1 - y0);
	float zd = (z - (float)z0) / (z1 - z0);

	// z interpolation
	float i1 = values[x0][y0][z0] * (1 - zd) + values[x0][y0][z1] * zd;
	float i2 = values[x0][y1][z0] * (1 - zd) + values[x0][y1][z1] * zd;
	float j1 = values[x1][y0][z0] * (1 - zd) + values[x1][y0][z1] * zd;
	float j2 = values[x1][y1][z0] * (1 - zd) + values[x1][y1][z1] * zd;

	// along y...
	float w1 = i1 * (1 - yd) + i2 * yd;
	float w2 = j1 * (1 - yd) + j2 * yd;

	// and return after interpolation along x
	return w1 * (1 - xd) + w2 * xd;
}

void point_lattice::iterator_reset()
{
    iteratorx = iteratory = iteratorz = 0;
}

float point_lattice::iterator_get_next_value(unsigned long* current_x,
                                             unsigned long* current_y,
                                             unsigned long* current_z)
{
    if (!iterator_is_done()) {
        float current_value = values[iteratorx][iteratory][iteratorz];
        (*current_x) = iteratorx;// * voxel_xlen;
        (*current_y) = iteratory;// * voxel_ylen;
        (*current_z) = iteratorz;// * voxel_zlen;
        // Update the position of the pointers for the next iteration
        // We update the most nested index (z) for sure
        // and update y and x only if either if we have finished a row and/or
        // column
        iteratorz++;
        // if the row is finished then set it back to zero and update y too
        if (iteratorz >= nz) {
            iteratorz = 0;
            iteratory++;
            if(iteratory >= ny) {
                iteratory = 0;
                iteratorx++;
                // if it goes over the last time then !iterator_is_done()
                // will catch it on the next call
            }
        }
        return current_value;
    } else {
        std::cerr << "Invalid density lattice iterator access. "
            << "The iterator finished but next value was requested again\n";
        throw 1;
    }
}

bool point_lattice::iterator_is_done()
{
    // When we go over the outermost loop (i.e. x-axis)
    return (iteratorx >= nx);
}
