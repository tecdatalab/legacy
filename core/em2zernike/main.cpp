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

/*
    This is a demonstration of how the classes may be used to generate the
    3D Zernike descriptors from a given input binary file containing the 
    voxel grid representation of the object.

    Notice that in the below case, the file contains a cubic grid, i.e. the 
    x-, y-, and z-dimensions are equal (such setting should satisfy most needs).
*/

#include <cstring>
#include <cstdlib>
// ---- local includes ----
#include "ZernikeDescriptor.h"

// ---- std includes ----
#include <vector>
#include <complex>
#include <fstream>
#include <sstream>
#include <iostream>

#include <stdio.h>

extern "C"{
    #include "lib_vio.h"
}


int main (int argc, char** argv)
{

    if (argc != 5)
    {
        std::cout << "Usage: ZernikeMoments " <<
                     "mapfilename " <<
                     "<MaxOrder> " <<
                     "<Contour> " <<
                     "outfileprefix " << std::endl;
        return 0;
    }


    // .inv file name
    char  buf1[200];
    char  buf2[200];
    //char  buf3[200];
    char* ch_ptr = strrchr (argv[1], '.');
    char* invFName;
    char* momFName;
    //char* recFName;

    if (ch_ptr)
    {
        //strncpy (buf1, argv[1], ch_ptr - argv[1]); 
        //buf1[ch_ptr - argv[1]] = '\0';
        strcpy (buf1, argv[4]); 
        strcpy (buf2, argv[4]); 
        //buf2[ch_ptr - argv[1]] = '\0';
        //strncpy (buf3, argv[1], ch_ptr - argv[1]); 
        //buf3[ch_ptr - argv[1]] = '\0';
    }
    else 
    {
        fprintf (stderr, "No extension in input filename? : %s\n", argv[1]);
    }

    invFName = buf1;
    momFName = buf2;
    //recFName = buf3;
    
    // compute the zernike descriptors
    sscanf(argv[3], "%lf", &CONTOUR); //parse the contour
    ZernikeDescriptor<double, float> zd (argv[1], atoi (argv[2]));

    strcat (invFName, ".inv");   
    std::cout << "Saving invariants file: " << invFName << " \n";
    // save them into an .inv file
    zd.SaveInvariants (invFName);

    strcat (momFName, ".3dzm");   
    std::cout << "Saving moments file: " << momFName << " \n";
    zd.SaveZernikeMoments (momFName);

}
