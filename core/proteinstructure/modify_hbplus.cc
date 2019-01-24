/*
 * This program takes the output of the HBPLUS program, which is a "$$$$.h" file
 * The program then modifies the atomtypes to the standard CHARMM types
 * and outputs a .pdb.ms file
 * TO COMPILE: g++ -O2 -Wall modify_hbplus.cc -o FORMATHBPLUS
 */



#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

const char kBlankChars[] = " \t\n\r";


class atom
{
    public:
        int         anum;   // id of the atom
        string      atype;  // type of the atom
        int         rnum;   // residue id
        double      axyz[3]; // xyz coordinates
        string      residue; // name of residue
        string      chain;// chain id
        int         atmhet;

        atom(int m_num, string m_type, int m_rid, double* m_xyz, string m_residue,
            string m_chain, int m_atmhet)
        {
            anum = m_num;
            rnum = m_rid;
            atype = m_type;
            for(int i=0;i<3;i++)
            {
                axyz[i] = m_xyz[i];
            }
            residue = m_residue;
            chain = m_chain;
            atmhet = m_atmhet;
        }

        atom(){}
        ~atom(){}
};

/// Returns a string with leading/trailing characters of a set stripped
string trimmed(string const& str, char const* sepSet=kBlankChars)
{
    string::size_type const first = str.find_first_not_of(sepSet);
    return (first==string::npos )? string():str.substr(first, str.find_last_not_of(sepSet)-first+1);
}

string rtrimmed(string const& str, char const* sepSet)
{
    string::size_type const last = str.find_last_not_of(sepSet);
    return (last==string::npos)? string():str.substr(0, last+1);
}

string ltrimmed(string const& str, char const* sepSet)
{
    string::size_type const first = str.find_first_not_of(sepSet);
    return (first==string::npos)? string():str.substr(first);
}


/**
 * return the filename after removing the path
 *
 * @params path filepath
 * @return basename of the file
 */
string basename(const string& path)
{
    string::size_type idx = path.find_last_of("\\/");

    if (idx == path.length()-1)
    {
        // last character is a path separator; i.e. path is a directory
        cerr << "invalid input file name" << endl;
        exit(EXIT_FAILURE);
    }

    if (idx != string::npos)
    {
        return path.substr(idx+1);
    }
    else
    {
        return path;
    }
}

/**
 * return the filename after removing the path and extension
 *
 * @params path filepath
 * @return actual file name
 */
string get_file_name(const string& path)
{
    // remove directory path from the filename
    string in_file = basename(path);

    // TODO: check the filename ends pdb

    // strip ".h" from the filename
    string fname = in_file.erase(in_file.length()-2);

    return fname;
}

void read_protein(string infile, vector<atom>& vct_atom)
{
    ifstream ifs;
    ifs.open(infile.c_str());
    if(!ifs)
    {
        cerr << "No file" << endl;
        exit(EXIT_FAILURE);
    }

    string line;

    int atom_id, res_id;
    string res_name, atom_type;
    double d_xyz[3];
    string chain_id;

    int atm_het = 0;

    while(getline(ifs, line))
    {
        if((line.substr(0, 4)).compare("ATOM") == 0 || (line.substr(0, 6)).compare("HETATM") == 0)
        {
            if((line.substr(0, 4)).compare("ATOM") == 0)
                atm_het = 1;
            else
                atm_het = 2;
                
            atom_id = atoi((line.substr(6, 5)).c_str());
            atom_type = trimmed(line.substr(12, 5), kBlankChars); // hydrogen atoms have 4 letter identification
            res_name = trimmed(line.substr(17, 3), kBlankChars);
            res_id = atoi((line.substr(22, 4)).c_str());
            chain_id = trimmed(line.substr(21, 1), kBlankChars);
            
            d_xyz[0] = atof((line.substr(30, 8)).c_str());
            d_xyz[1] = atof((line.substr(38, 8)).c_str());
            d_xyz[2] = atof((line.substr(46, 8)).c_str());

		/*
		 * This piece of code was originally in Vish's program, however it's effectively removing these
		 * atoms from the output
		 * This has been commented out unless we find that it's really necessary to remove them
            if(atom_type == "HE21" || atom_type == "HE22" || atom_type == "HD21" || atom_type == "HD22" ||
                atom_type == "H" || atom_type == "HZ1" || atom_type == "HZ2" || atom_type == "HZ3" ||
                atom_type == "HH21" || atom_type == "HH22" || atom_type == "HH11" || atom_type == "HH12")
                continue;
		*/

            if(atom_type == "1H" || atom_type == "2H" || atom_type == "3H")
                atom_type = "H";
            if(atom_type == "HE2")
                atom_type = "H"; // HIS
            if(atom_type == "HG")
                atom_type = "H";// CYS
            if(atom_type == "1HZ")
                atom_type = "HZ1";
            if(atom_type == "2HZ")
                atom_type = "HZ2";
            if(atom_type == "3HZ")
                atom_type = "HZ3";
            if(atom_type == "1HD2")
                atom_type = "HD21";
            if(atom_type == "2HD2")
                atom_type = "HD22";
            if(atom_type == "1HE2")
                atom_type = "HE21";
            if(atom_type == "2HE2")
                atom_type = "HE22";
            if(atom_type == "1HH1")
                atom_type = "HH11";
            if(atom_type == "2HH1")
                atom_type = "HH12";
            if(atom_type == "1HH2")
                atom_type = "HH21";
            if(atom_type == "2HH2")
                atom_type = "HH22";

            // create atom object
            atom n_atom(atom_id, atom_type, res_id, d_xyz, res_name, chain_id, atm_het);

            vct_atom.push_back(n_atom);
        }
    }// end while

    ifs.close();
}


void write_pdb(string fname, vector<atom>& P)
{
    string pdb_file_name = fname;
    FILE *pdb_file ;
    if((pdb_file = fopen( pdb_file_name.c_str(), "w")) == NULL)
    {
        printf( "This file could not be opened.\nDying\n\n" ) ;
        exit(EXIT_FAILURE);
    }
    cerr << "Writing data to " << pdb_file_name << endl;

    //fprintf( pdb_file, "ATOM  %5d %4s %3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f              %1s %2d\n",
    for(size_t i=0;i<P.size();i++)
    {
        if(P[i].atmhet == 1)
	{
		// atom types occupy 4 characters, however, when they have 3 characters or less the first one is used
		// as a space. The first character positions is only used if we actually use all 4 characters.
		// Given this fact we need to have two printing positions for the atom type
	    const char* print_template = P[i].atype.length() == 4 ?
	    	"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f\n":
	    	"ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f\n";
            fprintf(pdb_file, print_template,
                P[i].anum, P[i].atype.c_str(), P[i].residue.c_str(), P[i].chain.c_str(), P[i].rnum, P[i].axyz[0], P[i].axyz[1], P[i].axyz[2]);
	}
    }
    fclose(pdb_file);
}


int main(int argc, char** argv)
{
    if(argv[1] == NULL || argv[2] == NULL)
    {
        cerr << "Program pdb_file outfile" << endl;
        exit(EXIT_FAILURE);
    }

    string infile(argv[1]);
    string outfile(argv[2]);

    vector<atom> P;
    read_protein(infile, P);

    write_pdb(outfile, P);

    exit(EXIT_SUCCESS);
}

