#ifndef _ZDOCK_H_
#define _ZDOCK_H_

using namespace std;

#include "atom.h"
#include <vector>

class zdock
{
	public:
		int TVECT[3];
		double CANGLE[3];
		double SCORE;

		zdock(double* t1, int* t2, double score)
		{
			SCORE = score;
			for(int i=0;i<3;i++)
			{
				CANGLE[i] = t1[i];
			}
			for(int i=0;i<3;i++)
			{
				TVECT[i] = t2[i];
			}
		}
		
		zdock(){}
		~zdock(){}
};

class zdata
{
	public:
		double RANDVECT[3];
		double RVECT[3];
		double LVECT[3];
		double SPACING;
		int N;
		vector<zdock> ZD;
		zdata(){}
		~zdata(){}
};

#endif

