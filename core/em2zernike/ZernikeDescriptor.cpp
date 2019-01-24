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

extern "C"{
    #include "lib_vio.h"
}

template<class T, class TIn>
ZernikeDescriptor<T, TIn>::ZernikeDescriptor (const char* _rawName, int _order) : order_ (_order)
{
    printf("reading gridfile...\n");
    //voxels_ =  ReadGrid (_rawName, dim_);
    voxels_ =  ReadMap (_rawName, dim_);
    printf("computing...\n");

    // scale + translation normalization
    ComputeNormalization ();
    puts("finished ComputeNormalization");
    NormalizeGrid ();
    puts("finished NormalizeGrid");

    ComputeMoments ();
    puts("finished ComputeMoments");
    ComputeInvariants ();
    puts("finished ComputeInvariants");

//#ifdef	DEBUG_ZD	
		cout 	<< "\nxCOG: " << xCOG_
					<< "\nyCOG: " << yCOG_ 
					<< "\nzCOG: " << zCOG_ 
					<< "\nscale: "<< scale_ << endl;
//#endif
		
}



template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeMoments ()
{
    puts("initing GM");
    gm_.Init (voxels_, dim_, dim_, dim_, xCOG_, yCOG_, zCOG_, scale_, order_);

    // Zernike moments
    puts("initing ZM");
    zm_.Init (order_, gm_);
    puts("computing ZM");
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
    T sqrLen;

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

                    sqrLen = point[0]*point[0] + point[1]*point[1] + point[2]*point[2];
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
    puts("Constructing geometric moments...");
    ScaledGeometricalMoments<T, T> gm (voxels_, dim_, 0.0, 0.0, 0.0, 1.0);
    puts("   done.");
    
    // compute the geometrical transform for no translation and scaling, first
    // to get the 0'th and 1'st order properties of the function
    //gm.Compute (); 

    // 0'th order moments -> normalization
    // 1'st order moments -> center of gravity
    zeroMoment_ = gm.GetMoment (0, 0, 0);
    xCOG_ = gm.GetMoment (1, 0, 0) / zeroMoment_;
    yCOG_ = gm.GetMoment (0, 1, 0) / zeroMoment_;
    zCOG_ = gm.GetMoment (0, 0, 1) / zeroMoment_;
    puts("got 0/1th moments");
    // scaling, so that the function gets mapped into the unit sphere
    
    T recScale = 2.0 * ComputeScale_RadiusVar (voxels_, dim_, xCOG_, yCOG_, zCOG_);
    puts("got raduisvar");
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
    //char temp;

    // read electrostatic potential or hydrophobic grid values
    T tmp;

    size_t numpts = 0;
    //int nonzero = 0;

    while (infile >> tmp) {

        if(numpts%(128*(2<<20)) == 0)
        {
            tempGrid.reserve(numpts+128*(2<<20));
        }

        tempGrid.push_back(tmp);
        numpts++;
    }
    

    //int d = tempGrid.size ();
    size_t d = numpts;
    double f = pow ((double)d, 1.0/3.0);
    _dim_ = (int)floor (f+0.5);
    
    printf("File: %lu Grid: %d\n", d, _dim_);

    /*if (d != (size_t)(_dim_*_dim_*_dim_))
        puts("d != dim^3");*/

    d = _dim_*_dim_*_dim_;
        
    T* result = new T [d];
    for (size_t i=0; i<d; i++)
    {
        result[i] = 1000*tempGrid[i];
    }

    //size_t xdim=_dim_, ydim=_dim_/*, zdim=_dim_*/; //can be autooptimized away...
/*
    SparseGridT sgt;
    //int cc = 0;
    for(size_t i = 0; i < d; i++)
    {
        if(tempGrid[i] != 0)
        {
            GridPoint gp;
            //turn the index into grid coords
            gp.x = ((i%(xdim*ydim))%xdim);
            gp.y = ((i%(xdim*ydim))/xdim);
            gp.z = (i/(xdim*ydim));
            //fprintf(stderr, "point: %lu %lu %lu\n",gp.x, gp.y, gp.z);
            sgt[gp] = tempGrid[i];
            //if(++cc > 100)
            //    break;
        }
    }

    sparsegrid_ = sgt;
*/
    return result;
}

/**
 * Reads the grid from a binary file containing the float grid values. The
 * TIn type tells waht is the precision of the grid. In this implementstion
 * it is assumed that the dimensions of the grid are equal along each axis.
 * 
 * THIS ONE IS WITH MAP FILES
 */

template<class T, class TIn>
T* ZernikeDescriptor<T, TIn>::ReadMap (const char* _fname, int& _dim_)
{
    double porigx, porigy, porigz;  
    double *phi;
    double *pphi;
    unsigned pextx, pexty, pextz;
    unsigned long nvox;
    int ordermode = 7, swapmode, cubic = 1, orom = 1;
    double pwidth, widthx, widthy, widthz;
    double alpha, beta, gamma;
    int nxstart=0, nystart=0, nzstart=0;
    int ncstart=0, nrstart=0, nsstart=0;
    unsigned nc, nr, ns;
    unsigned nx, ny, nz;
    int currext;
    double xorigin, yorigin, zorigin;
    char ac='X', ar='Y', as='Z';

    size_t numpts = 0;

    char fname_nonconst[200];
    strcpy(fname_nonconst, _fname);

    /************** process input map ******************/ 
  
    
    /* MRC or CCP4 format with automatically filled cubic grid parameters */
    
    read_mrc(fname_nonconst, &orom, &cubic, &ordermode, &nc, &nr, &ns, &ncstart, &nrstart, &nsstart, &widthx, &widthy, &widthz, &xorigin, &yorigin, &zorigin, &alpha, &beta, &gamma, &phi);
    permute_map(ordermode, nc, nr, ns, &nx, &ny, &nz, ncstart, nrstart, nsstart, &nxstart, &nystart, &nzstart, phi, &pphi);
    permute_print_info(ordermode,nc,nr,ns,nx,ny,nz,ncstart,nrstart,nsstart,nxstart,nystart,nzstart);
    assert_cubic_map(orom, cubic, alpha, beta, gamma, widthx, widthy, widthz, nx, ny, nz, nxstart, nystart, nzstart, xorigin, yorigin, zorigin, &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);

    numpts = pextx < pexty ? pexty : pextx;
    numpts = pextz < numpts ? numpts : pextz;
    _dim_ = numpts;
    numpts = numpts*numpts*numpts;

    //int d = tempGrid.size ();
    size_t d = numpts;
    //double f = pow ((double)d, 1.0/3.0);
    //_dim_ = (int)floor (f+0.5);
    
    printf("File: %lu Grid: %d\n", d, _dim_);

    T* result = new T[d];

    size_t idx = 0, k = 0;

    //printf("PET: (%d, %d, %d)\n", pextx, pexty, pextz);

    /*FILE* ff = fopen("emd_5716.grid","r");
    FILE* gg = fopen("emd_5716.situs","r");
    int REF = -1, XX = -1;
    double REG = -1.0;
    double PP;
    fscanf(gg, "%lf", &REG);
    fscanf(gg, "%lf", &REG);
    fscanf(gg, "%lf", &REG);
    fscanf(gg, "%lf", &REG);
    fscanf(gg, "%lf", &REG);
    fscanf(gg, "%lf", &REG);
    fscanf(gg, "%lf", &REG);*/

    for(size_t z = 0; z < _dim_; z++)
        for(size_t y = 0; y < _dim_; y++)
            for(size_t x = 0; x < _dim_; x++)
            {
                if((x>=pextx) || (y>=pexty) || (z>=pextz))
                {
                    result[idx++] = 0;
                    //XX=0;
                    //PP=-1.0;
                    //abort();
                }
                else if(pphi[k]<=CONTOUR)
                {
                    //PP = pphi[k];
                    result[idx++] = 0;
                    k++;
                   // XX=1;
                }
                else
                {
                    //PP = pphi[k];
                    result[idx++] = /*(pphi[k++]**/1000;
                    //printf(" %f \n", pphi[k-1]);
                    k++;
                    //XX=1;
                }
               // fprintf(ff,"%d ", (result[idx-1]>0) ? 1 : 0);
               /* fscanf(ff,"%d",&REF);
                fscanf(gg,"%lf",&REG);
                if(REF !=(result[idx-1]>0))
                {
                    printf("INCONSISTENCY: XX:%d IDX:%lu K:%lu VAL:%lf REF:%d US:%d COORD:(%lu, %lu, %lu) SITUS:%lf PP:%lf CONTOUR:%lf\n", XX, idx-1, k-1, pphi[k-1], REF, (int)result[idx-1], x, y, z, REG, PP, CONTOUR);
                    abort();
                }*/
            }
            //fclose(ff);
            //fclose(gg);
    //printf("%lu, %lu\n", k, idx);
    //free(pphi);
    //free(phi);

    /************** write output map ******************/ 
    /* Situs or MRC/CCP4 format output determined by file extension */
    //write_vol(argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);

    return result;
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
        //int l0 = n % 2;
        int li = 0;

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
void ZernikeDescriptor<T,TIn>::SaveZernikeMoments (const char* _fName)
{
    std::ofstream outfile (_fName, std::ios_base::binary | std::ios_base::out);

    for (int n=0; n<order_+1; ++n)
    {
        for (int l = n % 2; l<=n; l+=2)
        {
            for (int m=-l; m<=l; ++m)
            {
                ComplexT moment = zm_.GetMoment (n, l, m);
                outfile << n << ' ' << l << ' ' << m << ' ' << moment << std::endl;
            }
        }
    }
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
