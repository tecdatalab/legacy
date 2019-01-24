#include "read.h"


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

    // strip ".pdb"/".pqr" from the filename
    string fname = in_file.erase(in_file.length()-4);

    return fname;
}


/**
 * Read PDB informtion from file
 * Parse only "ATOM" information.
 *
 * @see    process_files
 * @params infile name of pdb file
 * @params pdbdata stores pdb coordinates and atom types
 * @return void
 */

void read_protein(string infile, pdb& pdbdata)
{
    ifstream ifs;
    ifs.open(infile.c_str());
    if(!ifs)
    {
        cerr << "Could not open file: " << infile << endl;
        exit(EXIT_FAILURE);
    }

    string line;

    int atom_id, res_id;
    string chain_id;
    string res_name, atom_type, achar = "";
    double d_xyz[3];


    vector<atom> vct_atom;

    while(getline(ifs, line))
    {
        if((line.substr(0, 4)).compare("ATOM") == 0)
        {
            chain_id = trimmed(line.substr(21, 1), kBlankChars);
    
            atom_id = atoi((line.substr(6, 5)).c_str());
            atom_type = trimmed(line.substr(12, 4), kBlankChars); // hydrogen atoms have 4 letter identification
            res_name = trimmed(line.substr(17, 3), kBlankChars);
            
            res_id = atoi((line.substr(22, 4)).c_str());
            //cout << res_id << endl;
            d_xyz[0] = atof((line.substr(30, 8)).c_str());
            d_xyz[1] = atof((line.substr(38, 8)).c_str());
            d_xyz[2] = atof((line.substr(46, 8)).c_str());


            // create atom object
            atom n_atom(atom_id, atom_type, res_id, d_xyz, res_name, achar, chain_id);

            vct_atom.push_back(n_atom);
        }
    }// end while


    string fname = get_file_name(infile);
    
    pdbdata.atoms = vct_atom;
    pdbdata.pname = fname;

    ifs.close();
}



void read_visgrid_output(string infile, vector<vector<int> >& V)
{
    ifstream ifs;
    ifs.open(infile.c_str());
    if(!ifs)
    {
        cerr << "Could not open file: " << infile << endl;
        exit(EXIT_FAILURE);
    }


    string line;
    int d;
    
    while(getline(ifs, line))
    {
        vector<int> X;
        istringstream iss(line);
        while(iss >> d)
        {
            X.push_back(d);
        }
        V.push_back(X);
    }
    ifs.close();

}


