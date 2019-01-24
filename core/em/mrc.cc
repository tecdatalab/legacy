#include "mrc.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

mrc_reader::mrc_reader(string filename)
{
	mrc_file.open(filename.c_str(), ios::in | ios::binary);
    if(!mrc_file.is_open()) {
        // just a generic exception number
        cerr << "Could not open MRC file" << endl;
        throw 1;
    }
    test_and_set_endianness();
    // Read the header parts that we're interested in
    // Total columns(x), rows(y) and sections (z)
    this->nx = read_header_word();
    this->ny = read_header_word();
    this->nz = read_header_word();
    // Integral type used to store densities
    this->mode = read_header_word();
    if (mode > 2) {
      // We only handle modes 0, 1 and 2
      cerr << "Only map modes 0, 1 and 2 are supported. Mode "
           << mode << " was found" << endl;
      throw 1;
    }
    // Start x,y,z positions
    this->nxstart = read_header_word();
    this->nystart = read_header_word();
    this->nzstart = read_header_word();
    // Grid size
    this->mx = read_header_word();
    this->my = read_header_word();
    this->mz = read_header_word();
    // Cell dimensions
    this->xlen = read_four_byte_float();
    this->ylen = read_four_byte_float();
    this->zlen = read_four_byte_float();
    // Ignore 3 cell angle words and the 3 words corresponding to axes
    // which means that the pointer should be moved right before word #20
    mrc_file.seekg(19 * sizeof(mrcword_t));
    // Densitiy ranges
    this->dmin = read_four_byte_float();
    this->dmax = read_four_byte_float();
    this->dmean = read_four_byte_float();
    // Jump before the origin (word #49)
    mrc_file.seekg(49 * sizeof(mrcword_t));
    // Origin (transformations)
    this->xorigin = read_four_byte_float();
    this->yorigin = read_four_byte_float();
    this->zorigin = read_four_byte_float();

// Make a full copy of the header in case it needs to be rewritten.
    mrc_file.seekg(0);
    mrc_file.read((char*)header, 256 * sizeof(mrcword_t));

/*    
    cout << "Reading MRC file " << filename << "..." << endl;
    cout << "NX=" << nx << " NY=" << ny << " NZ=" << nz
        << " Mode=" << mode << endl;
    cout << "NXSTART=" << nxstart << " NYSTART=" << nystart
        << " NZSTART=" << nzstart << endl;
    cout << "MX=" << mx << " MY=" << my << " MZ=" << mz << endl;
    cout << "XLEN=" << xlen << " YLEN=" << ylen << " ZLEN=" << zlen << endl;
    cout << "DMIN=" << dmin << " DMAX=" << dmax
        << " DMEAN=" << dmean << endl;

    cout << "XORIGIN=" << xorigin << " YORIGIN=" << yorigin
        << " ZORIGIN=" << zorigin << endl;
    cout << "Voxel size x: " << (xlen / nx) << " y: " << (ylen / ny)
        << " z: " << (zlen / nz) << endl;
*/    
    // Move to the end of the header
    mrc_file.seekg(MRC_HEADER_SIZE);
/*    
    long density_start = mrc_file.tellg();
    mrc_file.seekg (0, ios::end);
    long density_end = mrc_file.tellg();
    
    cout << "Total voxels: " << ((density_end - density_start) / 4) << endl;
    cout << "Expected total voxels: " << (nx * ny * nz) << endl;
*/
}

mrc_reader::~mrc_reader()
{
    mrc_file.close();
}

void mrc_reader::test_and_set_endianness()
{
  // Make sure we're at the beginning to read the NX
  mrc_file.seekg(0);
  // The first header word is the number of columns.
  // We assume that the lowest number is the correct read direction
  mrcword_t regular_nx = read_header_word(false);
  mrcword_t reversed_nx = read_header_word(true);
  this->is_endianness_reversed = reversed_nx < regular_nx;
  // Rewind back to the beginning
  mrc_file.seekg(0);
}

mrcword_t mrc_reader::read_header_word()
{
  return read_header_word(is_endianness_reversed);
}

mrcword_t mrc_reader::read_header_word(bool is_reversed)
{
  mrcword_t value_read;
  mrc_file.read((char*)&value_read, sizeof(mrcword_t));
  if (is_reversed) {
    swap_four_bytes((unsigned char*)&value_read);
  }
  return value_read;
}

float mrc_reader::read_four_byte_float()
{
  float value_read;
  mrc_file.read((char*)&value_read, sizeof(mrcword_t));
  if (is_endianness_reversed) {
    swap_four_bytes((unsigned char*)&value_read);
  }
  return value_read;
}

float mrc_reader::read_mode0_density()
{
  char value_read;
  mrc_file.read((char*)&value_read, sizeof(char));
  return (float)value_read;
}

float mrc_reader::read_mode1_density()
{
  short value_read;
  mrc_file.read((char*)&value_read, sizeof(short));
  if (is_endianness_reversed) {
    swap_two_bytes((unsigned char*)&value_read);
  }
  return (float)value_read;
}

float mrc_reader::read_mode2_density()
{
  // It's the same as the header floats but for readability
  // we create a function just for this purpose.
  return read_four_byte_float();
}

void mrc_reader::swap_four_bytes(unsigned char* value)
{
  unsigned char* swap_pointer = value;
  unsigned char swap_value = swap_pointer[0];
  swap_pointer[0] = swap_pointer[3];
  swap_pointer[3] = swap_value;
  swap_value = swap_pointer[1];
  swap_pointer[1] = swap_pointer[2];
  swap_pointer[2] = swap_value;
}

void mrc_reader::swap_two_bytes(unsigned char* value)
{
  unsigned char* swap_pointer = value;
  unsigned char swap_value = swap_pointer[0];
  swap_pointer[0] = swap_pointer[1];
  swap_pointer[1] = swap_value;
}
density_lattice* mrc_reader::get_densities()
{
    density_lattice* densities =
        new density_lattice(nx, ny, nz, xorigin, yorigin, zorigin,
                            (xlen / nx), (ylen / ny), (zlen / nz));
    // Move to the end of the header
    mrc_file.seekg(MRC_HEADER_SIZE);
    
    for (mrcword_t z = 0; z < nz; z++) {
        for (mrcword_t y = 0; y < ny; y++) {
            for (mrcword_t x = 0; x < nx; x++) {
              float current_density;
              if (this->mode == 0) {
                current_density = read_mode0_density();
              } else if (this->mode == 1) {
                current_density = read_mode1_density();
              } else { // mode 2 (an exception raised earlier otherwise)
                current_density = read_mode2_density();
              }
              densities->set(current_density, x, y, z);
            }
        }
    }
    return densities;
}

void mrc_reader::set_origin(float new_x, float new_y, float new_z)
{
  xorigin = new_x;
  yorigin = new_y;
  zorigin = new_z;
  // Adjust the raw header too, in case it's written back to a file.
  // Floats are 4 bytes, and that's the size expected from the MRC spec.
  memcpy(header + 49, &new_x, sizeof(mrcword_t));
  memcpy(header + 50, &new_y, sizeof(mrcword_t));
  memcpy(header + 51, &new_z, sizeof(mrcword_t));
}

void mrc_reader::write(string filename)
{
  ofstream output(filename.c_str(), ios::out | ios::binary | ios::trunc);
  if(output.fail()) {
    cerr << "Could not write to MRC file: " << filename << endl;
  } else {
    output.write((char*)header, 256 * sizeof(mrcword_t));
    mrc_file.seekg(MRC_HEADER_SIZE);
    // Not the most efficient, but easy to read to avoid buffering
    mrcword_t total_voxels = nx * ny * nz;
    float current_density;
    for(mrcword_t voxel = 0; voxel < total_voxels; voxel++) {
      // TODO; This assumes Mode 2 representation. Change it to handle
      // Modes 0 and 1 (the same as the read operation does now).
      mrc_file.read((char*)&current_density, sizeof(float));
      output.write((char*)&current_density, sizeof(float));
    }
    output.close();
    if(output.fail()) {
      cerr << "Could not write to MRC file: " << filename << endl;
    }
  }
}
