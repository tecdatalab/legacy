#include <iostream>
#include "pdb.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		cout << "Usage: reader <pdbFile>" << endl;
		return 1;
	}
	else
	{
		pdb pdbdata;
		read_protein(argv[1], pdbdata);
		cout << pdbdata.pname << " Atoms " << pdbdata.atoms.size() << endl;
		return 0;
	}
}
