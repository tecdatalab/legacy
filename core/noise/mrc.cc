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
	// Read the header parts that we're interested in
	// Total columns(x), rows(y) and sections (z)
	mrc_file.read((char*)&nx, sizeof(mrcword_t));
	mrc_file.read((char*)&ny, sizeof(mrcword_t));
	mrc_file.read((char*)&nz, sizeof(mrcword_t));
	// Integral type used to store densities
	mrc_file.read((char*)&mode, sizeof(mrcword_t));
	// Start x,y,z positions
	mrc_file.read((char*)&nxstart, sizeof(mrcword_t));
	mrc_file.read((char*)&nystart, sizeof(mrcword_t));
	mrc_file.read((char*)&nzstart, sizeof(mrcword_t));
	// Grid size
	mrc_file.read((char*)&mx, sizeof(mrcword_t));
	mrc_file.read((char*)&my, sizeof(mrcword_t));
	mrc_file.read((char*)&mz, sizeof(mrcword_t));
	// Cell dimensions
	mrc_file.read((char*)&xlen, sizeof(mrcword_t));
	mrc_file.read((char*)&ylen, sizeof(mrcword_t));
	mrc_file.read((char*)&zlen, sizeof(mrcword_t));
	// Ignore 3 cell angle words and the 3 words corresponding to axes
	// which means that the pointer should be moved right before word #20
	mrc_file.seekg(19 * sizeof(mrcword_t));
	// Densitiy ranges
	mrc_file.read((char*)&dmin, sizeof(mrcword_t));
	mrc_file.read((char*)&dmax, sizeof(mrcword_t));
	mrc_file.read((char*)&dmean, sizeof(mrcword_t));
	// Jump before the origin (word #49)
	mrc_file.seekg(49 * sizeof(mrcword_t));
	// Origin (transformations)
	mrc_file.read((char*)&xorigin, sizeof(mrcword_t));
	mrc_file.read((char*)&yorigin, sizeof(mrcword_t));
	mrc_file.read((char*)&zorigin, sizeof(mrcword_t));

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

density_lattice* mrc_reader::get_densities()
{
	density_lattice* densities =
			new density_lattice(nx, ny, nz, xorigin, yorigin, zorigin,
					(xlen / nx), (ylen / ny), (zlen / nz));
	// Move to the end of the header
	mrc_file.seekg(MRC_HEADER_SIZE);

	float current_density;
	for (mrcword_t z = 0; z < nz; z++) {
		for (mrcword_t y = 0; y < ny; y++) {
			for (mrcword_t x = 0; x < nx; x++) {
				mrc_file.read((char*)&current_density, sizeof(mrcword_t));
				densities->set(current_density, x, y, z);
			}
		}
	}
	return densities;
}

/*
 * Overloaded Method
 * threshold > 0; should be the suggested contour level
 */
density_lattice* mrc_reader::get_densities(float threshold)
{
	density_lattice* densities =
			new density_lattice(nx, ny, nz, xorigin, yorigin, zorigin,
					(xlen / nx), (ylen / ny), (zlen / nz));
	// Move to the end of the header
	mrc_file.seekg(MRC_HEADER_SIZE);

	float current_density;
	for (mrcword_t z = 0; z < nz; z++) {
		for (mrcword_t y = 0; y < ny; y++) {
			for (mrcword_t x = 0; x < nx; x++) {
				mrc_file.read((char*)&current_density, sizeof(mrcword_t));

				//densities->set(current_density, x, y, z);

				//Updated to work with noise removal:
				if(current_density >= threshold) {
					densities->set(current_density, x, y, z);
				} else {
					densities->set(dmin, x, y, z);
				}
			}
		}
	}
	//Used for noise removal:
	densities->remove_noise(dmin);

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
			mrc_file.read((char*)&current_density, sizeof(float));
			output.write((char*)&current_density, sizeof(float));
		}
		output.close();
		if(output.fail()) {
			cerr << "Could not write to MRC file: " << filename << endl;
		}
	}
}

void mrc_reader::write(string filename, density_lattice* dl)
{
	ofstream output(filename.c_str(), ios::out | ios::binary | ios::trunc);
	if(output.fail()) {
		cerr << "Could not write to MRC file: " << filename << endl;
	} else {
		output.write((char*)header, 256 * sizeof(mrcword_t));
		//mrc_file.seekg(MRC_HEADER_SIZE);
		// Not the most efficient, but easy to read to avoid buffering
		//mrcword_t total_voxels = nx * ny * nz;
		float current_density;
		//for(mrcword_t voxel = 0; voxel < total_voxels; voxel++) {
		//	mrc_file.read((char*)&current_density, sizeof(float));
		for (mrcword_t z = 0; z < dl->size_z(); z++) {
				for (mrcword_t y = 0; y < dl->size_y(); y++) {
					for (mrcword_t x = 0; x < dl->size_x(); x++) {
						current_density = dl->get(x,y,z);
		                output.write((char*)&current_density, sizeof(float));
					}
				}
		}
		output.close();
		if(output.fail()) {
			cerr << "Could not write to MRC file: " << filename << endl;
		}
	}
}
