/* prototypes for read.cc */
#ifndef _READ_H_
#define _READ_H_

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

#include "pdb.h"
typedef double TRMTX[4][3];

const char kBlankChars[] = " \t\n\r";

void read_protein(string, pdb&);
void read_visgrid_output(string, vector<vector<int> >&);



#endif


