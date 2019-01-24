#ifndef _POINT_LATTICE_H_
#define _POINT_LATTICE_H_

/*
 * Base class for density_lattice
 * It abstracts what's related purely to the storage
 * and lenghts, and not specific things like correlation
 * computations.
 */
class point_lattice
{
  private:
    // Number of voxels along each axis
    unsigned long nx;
    unsigned long ny;
    unsigned long nz;
    // In real space, where the lattice start. For practical
    // purposes this is just a translation used when finding what voxel
    // needs to be read when a real coordinate is provided and the density
    // value needs to be found
    float xorigin;
    float yorigin;
    float zorigin;
    // Similar to the origin, required to find real-space coordinates.
    float voxel_xlen;
    float voxel_ylen;
    float voxel_zlen;
    float ***values;
    // Although a separate iterator class could be created,
    // for simplicity purposes the x,y,z of an iterator are kept.
    unsigned long iteratorx;
    unsigned long iteratory;
    unsigned long iteratorz;
  public:
    // Allocates the 3D lattice and initializes it to zeroes
    // Note that the lengths expected are per voxel, not total x, y, z
    // lengths like the ones read directly from the MRC file
    point_lattice(unsigned long in_nx, unsigned long in_ny,
                  unsigned long in_nz, float in_xorigin, float in_yorigin,
                  float in_zorigin, float in_voxel_xlen, float in_voxel_ylen,
                  float in_voxel_zlen);
    // Just release the lattice memory
    ~point_lattice();
	  // Accessors for voxel totals...
	  inline unsigned long size_x() {return nx;}
	  inline unsigned long size_y() {return ny;}
	  inline unsigned long size_z() {return nz;}
    // ...the origin
    inline float origin_x() {return xorigin;}
    inline float origin_y() {return yorigin;}
    inline float origin_z() {return zorigin;}
	  // ...and voxel sizes
	  inline float voxel_length_x() {return voxel_xlen;}
	  inline float voxel_length_y() {return voxel_ylen;}
	  inline float voxel_length_z() {return voxel_zlen;}
    // Returns density at voxel [x,y,z]
    inline float get(unsigned long x, unsigned long y, unsigned long z) {
      if (x < nx && y < ny && z < nz) {
        return values[x][y][z];
      } else {
        return 0.0;
      }
    }
    // Sets the density at voxel [x,y,z]
    inline void set(float value, unsigned long x, unsigned long y,
             unsigned long z) {
        values[x][y][z] = value;
    }

    // Using the origin and voxel length parameters, map a real coordinate
    // x, y, z to the voxel closest to it (after translating and scaling).
    // If the real coordinate falls outside the range of the lattice it will be
    // assigned a voxel closest to the edge of the grid.
    void real_coordinate_to_voxel(double real_x, double real_y, double real_z,
        unsigned long* voxel_x, unsigned long* voxel_y, unsigned long* voxel_z);
    // Symmetrical operation for voxel to real coordinate.
    void voxel_to_real_coordinate(unsigned long voxel_x, unsigned long voxel_y,
        unsigned long voxel_z, double* real_x, double* real_y, double* real_z);
    // Internally calls get_interpolated but it first scales down the real
    // coordinate and translates it to make it a coordinate in voxel space.
    float get_from_real_coordinate(double x, double y, double z);
	  // When the x,y,z coordinate looked for is not one of the exact
	  // lattice point, we need to do trilinear interpolation
	  // based on the values around the requested coordinates
	  float get_interpolated(float x, float y, float z);
	  // Bring all iterator indices to zero
    void iterator_reset();
    // Retrieve the current value and advance the iterator pointers
    float iterator_get_next_value(unsigned long* current_x,
                                  unsigned long* current_y,
                                  unsigned long* current_z);
	  // Alternate version that doesn't care about the coordinates
	  float iterator_get_next_value() {
		  unsigned long dummy_x, dummy_y, dummy_z;
		  return iterator_get_next_value(&dummy_x, &dummy_y, &dummy_z);
	  }
    // Check if we reached the end 
    bool iterator_is_done();
};

#endif
