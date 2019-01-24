#include <iostream>
#include "pdb.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		cout << "Usage: splitter <pdbFile>" << endl;
	}
	else
	{
		generate_chain_split_files(argv[1]);
	}
	return 0;
}
