#include <iostream>
#include <vector>
#include <cstdlib>
#include "pdb.h"
#include "atom.h"

using std::cerr;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char** argv)
{
	if(argv[1] == NULL)
	{
		cerr << "Usage: clashes file.pdb [out.pdb]" << endl;
		exit(EXIT_FAILURE);
	}

	string file(argv[1]);
	string outfile = argc > 2 ? string(argv[2]) : "";

	/* Read the chains */
	vector< vector<atom> > chains = read_chains(file);
	vector<atom> clashes = get_clashing_atoms(chains);
	vector< vector<atom> > towrite;
	towrite.push_back(clashes);

	write_complex(outfile, towrite);

	exit(EXIT_SUCCESS);
}
