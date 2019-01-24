/*--------------------------------------------------
 * main function to test reconstruction of zernike
 * discriptors. 
 * -------------------------------------------------*/

#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#include "ZernikeDescriptor.h"


template<class T, class TIn>
T* ReadGrid (const char* _fname, int& _dim_)
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
    //while (infile.read ((char*)(&temp), sizeof (TIn)))
    //{
    //    tempGrid.push_back ((T)temp);
    //}
   
	// read electrostatic potential or hydrophobic grid values
    T tmp;
    while (infile >> tmp) {
	tempGrid.push_back (tmp);
    }
    
    int d = tempGrid.size ();
    double f = pow ((double)d, 1.0/3.0);
    _dim_ = floor (f+0.5);
    
    printf("File: %d Grid: %d\n", d, _dim_);
    if (d != _dim_*_dim_*_dim_)
	    printf("David La???\n");
    d = _dim_*_dim_*_dim_;
	    
    T* result = new T [d];
    for (int i=0; i<d; ++i)
    {
    	    T f = tempGrid[i];
	    result[i] = 100*f;
    	    /*switch(n) {
    	    	case '1':       result[i] = 1000; break;
    	    	case '2':       result[i] = 666; break;
    	    	case '3':       result[i] = 333; break;
    	    	default: 	result[i] = 0;
    	    }*/
    }

    return result;
}

int main(int argc, char * argv[])
{

	cout << "in main\n";
	
	cout << "three arguments \n";
	ZernikeDescriptor<double, char> zd(argv[1], atoi( argv[2] ));
	
	char buf[100];
	char* ch_ptr = strrchr(argv[1], '.');
	char *invfile;
	if(ch_ptr)
	{
		strncpy (buf, argv[1], ch_ptr - argv[1]);
		buf[ch_ptr - argv[1]] = '\0';
	}
	
	invfile = buf;
	
	strcat(invfile, ".inv");
	
	zd.SaveInvariants(invfile);

	return 1;
}

