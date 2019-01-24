#ifndef _MRC_H_
#define _MRC_H_

/* MRC files are a standard format for electron density maps
 * The first 1024 bytes represent the header followed by the actual
 * densities.
 * 
 * 56 4-byte words compose the standard fields followed by
 * custom text labels that can be added (e.g. to express symmetry)
 *
 * This is a specification excerpt that describes those fields
 * (from http://ami.scripps.edu/software/mrctools/mrc_specification.php)
 * (Here's another description http://bio3d.colorado.edu/imod/doc/mrc_format.txt)
 * 1	NX       number of columns (fastest changing in map)
 * 2	NY       number of rows   
 * 3	NZ       number of sections (slowest changing in map)
 * 4	MODE     data type :
 *      0        image : signed 8-bit bytes range -128 to 127 (char)
 *      1        image : 16-bit halfwords (signed short)
 *      2        image : 32-bit reals (float)
 *      3        transform : complex 16-bit integers (2 * signed short)
 *      4        transform : complex 32-bit reals (2 * float)
 *      6        image : unsigned 16-bit range 0 to 65535 (unsigned short)
 * 5	NXSTART number of first column in map (Default = 0)
 * 6	NYSTART number of first row in map
 * 7	NZSTART number of first section in map
 * 8	MX       number of intervals along X
 * 9	MY       number of intervals along Y
 * 10	MZ       number of intervals along Z
 * 11-13	CELLA    cell dimensions in angstroms (xlen, ylen, zlen)
 * 14-16	CELLB    cell angles in degrees
 * 17	MAPC     axis corresp to cols (1,2,3 for X,Y,Z)
 * 18	MAPR     axis corresp to rows (1,2,3 for X,Y,Z)
 * 19	MAPS     axis corresp to sections (1,2,3 for X,Y,Z)
 * 20	DMIN     minimum density value
 * 21	DMAX     maximum density value
 * 22	DMEAN    mean density value
 * 23	ISPG     space group number 0 or 1 (default=0)
 * 24	NSYMBT   number of bytes used for symmetry data (0 or 80)
 * 25-48(49??)	EXTRA    extra space used for anything   - 0 by default
 * Note: There appears to be a discrepancy, it looks like the origin is 49-51
 * 50-52	ORIGIN   origin in X,Y,Z used for transforms
 * 53	MAP      character string 'MAP ' to identify file type
 * 54	MACHST   machine stamp
 * 55	RMS      rms deviation of map from mean density
 * 56	NLABL    number of labels being used
 * 57-256	LABEL(20,10) 10 80-character text labels
 *
 */

#include "density_lattice.h"

#include <fstream>

using namespace std;

#define MRC_HEADER_SIZE 1024
 
// The header contains a series of 4-byte values
typedef unsigned int mrcword_t;

//http://www.cplusplus.com/doc/tutorial/files/

class mrc_reader
{
 private:
  // File handle used to traverse the MRC files
  ifstream mrc_file;
  // Stores the raw header (values below are extracted from this array)
  mrcword_t header[256];
  // header values that we're interested in
  mrcword_t nx;
  mrcword_t ny;
  mrcword_t nz;
  mrcword_t mode;
  mrcword_t nxstart;
  mrcword_t nystart;
  mrcword_t nzstart;
  mrcword_t mx;
  mrcword_t my;
  mrcword_t mz;
  // These 3 lengths are total lengths so xlen / nx
  // is the cell dimension along the x-axis (similar for y and z)
  float xlen;
  float ylen;
  float zlen;
  float dmin;
  float dmax;
  float dmean;
  float xorigin;
  float yorigin;
  float zorigin;
  // By default we won't reverse words in the header but after calling
  // test_endianness_reversal
  bool is_endianness_reversed;
  // Set the endinanness flag by checking the first field of the header.
  // The first field represents the number of columns so it will interpret the
  // value using both endinanness standards and keep the one that yields the
  // lowest value. It's a heuristic but normally the number of columns is not
  // high, thus the incorrect endinanness should yield a pretty significant
  // difference (exposing the correct/incorrect).
  // Note: The file stream is reset to the beginning.
  void test_and_set_endianness();
  // Using the already open mrf_file input stream, and assuming the file pointer
  // is on the header section, read 4 bytes, taking into account the
  // is_endianness_reversed flag.
  mrcword_t read_header_word();
  // Overload that uses the given flag and not the instance variable
  mrcword_t read_header_word(bool is_reversed);
  float read_four_byte_float();
  // Auxiliary functions to swap independent of the original
  // data type (as long as the sizes are correct)
  void swap_four_bytes(unsigned char* value);
  void swap_two_bytes(unsigned char* value);
  // Since desnities have different sizes depending on the type, we
  // provide a function for each. All of them return a float generically,
  // since that's the representation used by other classes.
  float read_mode0_density();
  float read_mode1_density();
  float read_mode2_density();
 public:
  // Loads the header and leaves the file open for
	// the normal operation invoked subsequently: read_densities
	mrc_reader(string filename);
	// Closes the file handle opened by the constructor
	~mrc_reader();
  // Read the non-header part of the file and populate a
  // density lattice with the information
  density_lattice* get_densities();
  density_lattice* get_densities(float threshold);
  inline float get_x_origin() {return xorigin;}
  inline float get_y_origin() {return yorigin;}
  inline float get_z_origin() {return zorigin;}
  void set_origin(float new_x, float new_y, float new_z);
  // Writes the header in-memory (because it could have been modified)
  // but all the densities will be copied directly from the file read.
  void write(string filename); 
  void write(string filename, density_lattice* dl);
};

#endif
