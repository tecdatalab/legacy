/*
                                                                            
                          3D Zernike Moments                                
    Copyright (C) 2003 by Computer Graphics Group, University of Bonn       
           http://www.cg.cs.uni-bonn.de/project-pages/3dsearch/             
                                                                            
Code by Marcin Novotni:     marcin@cs.uni-bonn.de
       
for more information, see the paper:

@inproceedings{novotni-2003-3d,
    author = {M. Novotni and R. Klein},
    title = {3{D} {Z}ernike Descriptors for Content Based Shape Retrieval},
    booktitle = {The 8th ACM Symposium on Solid Modeling and Applications},
    pages = {216--225},
    year = {2003},
    month = {June},
    institution = {Universit\"{a}t Bonn},
    conference = {The 8th ACM Symposium on Solid Modeling and Applications, June 16-20, Seattle, WA}
}        
 *---------------------------------------------------------------------------* 
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/


template<class T, class TIn>
ZernikeDescriptor<T, TIn>::ZernikeDescriptor (T* _voxels, int _dim, int _order) :
    voxels_ (_voxels), dim_ (_dim), order_ (_order)
{
   
	  ComputeNormalization ();
    NormalizeGrid ();

    ComputeMoments ();
    ComputeInvariants ();
}

template<class T, class TIn>
ZernikeDescriptor<T, TIn>::ZernikeDescriptor (const char* _rawName, int _order) : order_ (_order)
{
    voxels_ =  ReadGrid (_rawName, dim_);

    // scale + translation normalization
    ComputeNormalization ();
    NormalizeGrid ();

    ComputeMoments ();
    ComputeInvariants ();

#ifdef	DEBUG_ZD	
		cout 	<< "\nxCOG: " << xCOG_
					<< "\nyCOG: " << yCOG_ 
					<< "\nzCOG: " << zCOG_ 
					<< "\nscale: "<< scale_ << endl;
#endif
		
}

/* 
 * modified Oct. 16, 2007 by Sael lee  
 * 
 * input center of gravity and radius interms of grid position
 * 
 * */
template<class T, class TIn>
ZernikeDescriptor<T, TIn>::ZernikeDescriptor (const char* _rawName, int _order, T _xCOG, T _yCOG, T _zCOG, T _scale) :
    order_ (_order), xCOG_ (_xCOG), yCOG_ (_yCOG), zCOG_ (_zCOG), scale_ (_scale)
{
	 	voxels_ =  ReadGrid (_rawName, dim_);
		NormalizeGrid ();

		ComputeMoments();
		ComputeInvariants ();

#ifdef	DEBUG_ZD	
		cout 	<< "\nxCOG: " << xCOG_
					<< "\nyCOG: " << yCOG_ 
					<< "\nzCOG: " << zCOG_ 
					<< "\nscale: "<< scale_ << endl;
#endif
	
}


template<class T, class TIn>
ZernikeDescriptor<T, TIn>::ZernikeDescriptor (T* _voxels, int _dim, int _order, T _xCOG, T _yCOG, T _zCOG, T _scale) :
    voxels_ (_voxels), dim_ (_dim), order_ (_order), xCOG_ (_xCOG), yCOG_ (_yCOG), zCOG_ (_zCOG), scale_ (_scale)
{
		
		NormalizeGrid ();

		ComputeMoments();
		ComputeInvariants ();

#ifdef	DEBUG_ZD	
		cout 	<< "\nxCOG: " << xCOG_
					<< "\nyCOG: " << yCOG_ 
					<< "\nzCOG: " << zCOG_ 
					<< "\nscale: "<< scale_ << endl;
#endif
	

}

template<class T, class TIn>
ZernikeDescriptor<T, TIn>::ZernikeDescriptor (int _order) :
    order_ (_order)
{
		zm_.Init(_order);
}

template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::Init(int _order)
{
		order_ = _order;
		zm_.Init(_order);
}

template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeLocal (const char* _rawName,T _xCOG, T _yCOG, T _zCOG, T _scale)
{
	
	 	voxels_ 	= ReadGrid (_rawName, dim_);
		xCOG_ 		= _xCOG;
	 	yCOG_ 		= _yCOG;
	 	zCOG_	 		= _zCOG;
	 	scale_		= _scale;

		NormalizeGrid ();

		ComputeMomentsLocal();
		ComputeInvariants ();

#ifdef	DEBUG_ZD	
		cout 	<< "\nxCOG: " << xCOG_
					<< "\nyCOG: " << yCOG_ 
					<< "\nzCOG: " << zCOG_ 
					<< "\nscale: "<< scale_ << endl;
#endif
	
}


template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeLocal (T* _voxels, int _dim,T _xCOG, T _yCOG, T _zCOG, T _scale)
{
	 	voxels_ 	= _voxels;
	 	dim_ 		= _dim;
	 	xCOG_ 		= _xCOG;
	 	yCOG_ 		= _yCOG;
	 	zCOG_	 	= _zCOG;
		scale_		= _scale;
	 
		NormalizeGrid ();
		
		ComputeMomentsLocal();
		ComputeInvariants ();

#ifdef	DEBUG_ZD	
		cout 	<< "\nxCOG: " << xCOG_
					<< "\nyCOG: " << yCOG_ 
					<< "\nzCOG: " << zCOG_ 
					<< "\nscale: "<< scale_ << endl;
#endif
	
}

template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeLocal (const char* _rawName)
{
	
	 	voxels_ 	= ReadGrid (_rawName, dim_);
		
		ComputeNormalization ();
	 	NormalizeGrid ();

		ComputeMomentsLocal();
		ComputeInvariants ();


#ifdef	DEBUG_ZD	
		cout 	<< "\nxCOG: " << xCOG_
					<< "\nyCOG: " << yCOG_ 
					<< "\nzCOG: " << zCOG_ 
					<< "\nscale: "<< scale_ << endl;
#endif
	
}


template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeLocal (T* _voxels, int _dim)
{
	 	voxels_ 	= _voxels;
	 	dim_ 		= _dim;
		
		ComputeNormalization ();
		NormalizeGrid ();
		
		ComputeMomentsLocal();
		ComputeInvariants ();	


#ifdef	DEBUG_ZD	
		cout 	<< "\nxCOG: " << xCOG_
					<< "\nyCOG: " << yCOG_ 
					<< "\nzCOG: " << zCOG_ 
					<< "\nscale: "<< scale_ << endl;
#endif
	
}


template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeMomentsLocal()
{	
	gm_.Init (voxels_, dim_, dim_, dim_, xCOG_, yCOG_, zCOG_, scale_, order_);
	zm_.Compute (gm_);
}

/*
 * END MODIFIED
 */


template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeMoments ()
{
    gm_.Init (voxels_, dim_, dim_, dim_, xCOG_, yCOG_, zCOG_, scale_, order_);
    //gm_.SetTransform (xCOG_, yCOG_, zCOG_, scale_);
    //gm_.Compute ();

    // Zernike moments
    zm_.Init (order_, gm_);
    zm_.Compute ();
}

/**
 * Cuts off the function : the object is mapped into the unit ball according to
 * the precomputed center of gravity and scaling factor. All the voxels remaining
 * outside the unit ball are set to zero.
 */
template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::NormalizeGrid ()
{
    T point[3];

    // it is easier to work with squared radius -> no sqrt required
    T radius = (T)1 / scale_;
    T sqrRadius = radius * radius;

    for (int x=0; x<dim_; ++x)
    {
        for (int y=0; y<dim_; ++y)
        {
            for (int z=0; z<dim_; ++z)
            {
                if (voxels_[(z*dim_ + y)*dim_ + x] != (T)0)
                {
                    point[0] = (T)x - xCOG_;
                    point[1] = (T)y - yCOG_;
                    point[2] = (T)z - zCOG_;

                    T sqrLen = point[0]*point[0] + point[1]*point[1] + point[2]*point[2];
                    if (sqrLen > sqrRadius)
                    {
                        voxels_[(z*dim_ + y)*dim_ + x] = 0.0;
                    }
                }
            }
        }
    }
}

/**
 * Center of gravity and a scaling factor is computed according to the geometrical
 * moments and a bounding sphere around the cog.
 */
template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeNormalization ()
{
    ScaledGeometricalMoments<T, T> gm (voxels_, dim_, 0.0, 0.0, 0.0, 1.0);
    
    // compute the geometrical transform for no translation and scaling, first
    // to get the 0'th and 1'st order properties of the function
    //gm.Compute (); 

    // 0'th order moments -> normalization
    // 1'st order moments -> center of gravity
    zeroMoment_ = gm.GetMoment (0, 0, 0);
    xCOG_ = gm.GetMoment (1, 0, 0) / zeroMoment_;
    yCOG_ = gm.GetMoment (0, 1, 0) / zeroMoment_;
    zCOG_ = gm.GetMoment (0, 0, 1) / zeroMoment_;

    // scaling, so that the function gets mapped into the unit sphere
    
    //T recScale = ComputeScale_BoundingSphere (voxels_, dim_, xCOG_, yCOG_, zCOG_);
    T recScale = 2.0 * ComputeScale_RadiusVar (voxels_, dim_, xCOG_, yCOG_, zCOG_);
    if (recScale == 0.0)
    {
        std::cerr << "\nNo voxels in grid!\n";
        exit (-1);
    }
    scale_ = (T)1 / recScale;                          
}

/**
 * Reads the grid from a binary file containing the float grid values. The
 * TIn type tells waht is the precision of the grid. In this implementstion
 * it is assumed that the dimensions of the grid are equal along each axis.
 */

template<class T, class TIn>
T* ZernikeDescriptor<T, TIn>::ReadGrid (const char* _fname, int& _dim_)
{
	using namespace std;
    ifstream infile (_fname,  std::ios_base::in);
    if (!infile.is_open())
    {
        std::cerr << "Cannot open " << _fname << " for reading.\n";
        exit (-1);
    }

    vector<T> tempGrid;
    char temp;


/*	
    // read the grid values
    while (infile.read ((char*)(&temp), sizeof (char)))
    {

        tempGrid.push_back ((T)(temp-48));
    }
  
*/	
	// read electrostatic potential or hydrophobic grid values
    T tmp;
    while (infile >> tmp) {
		tempGrid.push_back (tmp);
    }
    

    int d = tempGrid.size ();
    double f = pow ((double)d, 1.0/3.0);
    _dim_ = (int)floor (f+0.5);
    
    printf("File: %d Grid: %d\n", d, _dim_);
    if (d != _dim_*_dim_*_dim_)
	    printf("d != dim^3\n");
    d = _dim_*_dim_*_dim_;
	    
    T* result = new T [d];
    for (int i=0; i<d; i++)
    {
	    T f = tempGrid[i];
			result[i] = 1000*f;
    }

	/*
	//test
	int x =0;
	for(int i=0; i<_dim_; i++){
		for(int j=0; j<_dim_; j++){
			for(int k=0; k<_dim_; k++){
				cout << result[x++] << " ";
			}
			cout << endl;
		}
	}
	*/

    return result;
}



/*
template<class T, class TIn>
T* ZernikeDescriptor<T, TIn>::ReadGrid (const char* _fname, int& _dim_)
{
    std::ifstream infile (_fname, std::ios_base::binary | std::ios_base::in);
    if (!infile)
    {
        std::cerr << "Cannot open " << _fname << " for reading.\n";
        exit (-1);
    }

    vector<T> tempGrid;
    TIn temp;

    // read the grid values
    while (infile.read ((char*)(&temp), sizeof (TIn)))
    {
        tempGrid.push_back ((T)temp);
    }  

    int d = tempGrid.size ();
    double f = pow ((double)d, 1.0/3.0);
    _dim_ = floor (f+0.5);

    T* result = new T [d];
    for (int i=0; i<d; ++i)
    {
        result[i] = tempGrid[i];
    }

    return result;
}

*/


template<class T, class TIn>
double ZernikeDescriptor<T, TIn>::ComputeScale_BoundingSphere (T* _voxels, int _dim, T _xCOG, T _yCOG, T _zCOG)
{
    T max = (T)0;

    // the edge length of the voxel grid in voxel units
    int d = _dim;

    for (int x=0; x<d; ++x)
    {
        for (int y=0; y<d; ++y)
        {
            for (int z=0; z<d; ++z)
            {
                if (_voxels[(z + d * y) * d + x] != 0)
                {
                    T mx = (T)x - _xCOG;
                    T my = (T)y - _yCOG;
                    T mz = (T)z - _zCOG;
                    T temp = mx*mx + my*my + mz*mz;

                    if (temp > max)
                    {
                        max = temp;
                    }
                }
            }
        }
    }

    return std::sqrt (max);
}

template<class T, class TIn>
double ZernikeDescriptor<T, TIn>::ComputeScale_RadiusVar (T* _voxels, int _dim, T _xCOG, T _yCOG, T _zCOG)
{
    // the edge length of the voxel grid in voxel units
    int d = _dim;

    int nVoxels = 0;

    T sum = 0.0;

    for (int x=0; x<d; ++x)
    {
        for (int y=0; y<d; ++y)
        {
            for (int z=0; z<d; ++z)
            {
                if (_voxels[(z + d * y) * d + x] != 0)
                {
                    T mx = (T)x - _xCOG;
                    T my = (T)y - _yCOG;
                    T mz = (T)z - _zCOG;
                    T temp = mx*mx + my*my + mz*mz;

                    sum += temp;

                    nVoxels++;
                }
            }
        }
    }

    T retval = sqrt(sum/nVoxels);

    return retval;
}

// Added by DKlab
template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::Reconstruct (T3D& _grid, int _minN, int _maxN, int _minL, int _maxL)
{
// the scaling between the reconstruction and original grid
//    T fac = (T)(_grid.size ()) / (T)dim_;

	
		//initialize the globle variable for results	
		cout <<"initialize the  recon_result\n";
		_grid.resize(dim_);
		for(int i = 0; i < dim_; i++){
			_grid[i].resize(dim_);
			for(int j = 0; j < dim_; j++){
				_grid[i][j].resize(dim_);
			}
		}	

//		T fac = (T)(_grid.size ()) / (T)dim_;
   	T fac = 1.0; 
		zm_.Reconstruct (
				_grid,         // result grid
				xCOG_*fac,     // center of gravity properly scaled 
				yCOG_*fac, 
				zCOG_*fac, 
				scale_/fac,    // scaling factor 
				_minN, _maxN,  // min and max freq. components to be reconstructed
				_minL, _maxL); 

}



template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::Reconstruct (ComplexT3D& _grid, int _minN, int _maxN, int _minL, int _maxL)
{
    // the scaling between the reconstruction and original grid
    T fac = (T)(_grid.size ()) / (T)dim_;

    zm_.Reconstruct (_grid,         // result grid
                     xCOG_*fac,     // center of gravity properly scaled 
                     yCOG_*fac, 
                     zCOG_*fac, 
                     scale_/fac,    // scaling factor 
                     _minN, _maxN,  // min and max freq. components to be reconstructed
                     _minL, _maxL); 
}


/**
 * Computes the Zernike moment based invariants, i.e. the norms of vectors with
 * components of Z_nl^m with m being the running index.
 */
template<class T, class TIn>
void ZernikeDescriptor<T,TIn>::ComputeInvariants ()
{
    //invariants_.resize (order_ + 1);
    invariants_.clear ();
    for (int n=0; n<order_+1; ++n)
    {
        //invariants_[n].resize (n/2 + 1);

        T sum = (T)0;
        int l0 = n % 2, li = 0;

        for (int l = n % 2; l<=n; ++li, l+=2)
        {
            for (int m=-l; m<=l; ++m)
            {
                ComplexT moment = zm_.GetMoment (n, l, m);
                sum += std::norm (moment);
            }

            invariants_.push_back (sqrt (sum));
            //invariants_[n][li] = std::sqrt (sum);
        }
    }
}

template<class T, class TIn>
void ZernikeDescriptor<T,TIn>::SaveInvariants (vector<float>&  inv)
{

    int dim = invariants_.size ();
		inv.resize(dim+1);
		
    inv[0] = (float)dim;
    for (int i=0; i<dim; ++i)
    {
        inv[i+1] = invariants_[i];
    }

		//SaveInvariants("testinzd.inv");
}

template<class T, class TIn>
void ZernikeDescriptor<T,TIn>::SaveInvariants (const char* _fName)
{
    std::ofstream outfile (_fName, std::ios_base::binary | std::ios_base::out);

    float temp;

    int dim = invariants_.size ();
    //outfile.write ((char*)(&dim), sizeof(int));
    outfile << dim << std::endl;
    for (int i=0; i<dim; ++i)
    {
        temp = invariants_[i];
	temp /= 10;
        //outfile.write ((char*)(&temp), sizeof(float));
	outfile << temp << std::endl;
    }
}

/*
template<class T, class TIn>
void ZernikeDescriptor<T,TIn>::SaveInvariants (const char* _fName)
{
    std::ofstream outfile (_fName, std::ios_base::binary | std::ios_base::out);

    float temp;

    int dim = invariants_.size ();
    outfile.write ((char*)(&dim), sizeof(int));
    
    for (int i=0; i<dim; ++i)
    {
        temp = invariants_[i];
        outfile.write ((char*)(&temp), sizeof(float));    
    }
}

*/




