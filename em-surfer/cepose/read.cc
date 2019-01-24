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

void read_protein(string infile, string ch, pdb& pdbdata)
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
            if(ch != "_")
            {
                if(chain_id != ch)
                    continue;
            }
    
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


void get_filename_chain(string line, string& protein, string& chain)
{
    istringstream iss(line);
    string s1, s2, s3;
    iss >> s1 >> s2 >> s3;

    size_t pos;
    pos = s3.find(":");

    protein = s3.substr(0, pos);
    chain = s3.substr(pos+1);
}

//X2 = ( 0.411764)*X1 + ( 0.442401)*Y1 + (-0.796700)*Z1 + (   62.513756)
void get_coords(string line, double* X)
{
    X[0] = atof(line.substr(6, 9).c_str());
    X[1] = atof(line.substr(23, 9).c_str());
    X[2] = atof(line.substr(40, 9).c_str());
    X[3] = atof(line.substr(57, 12).c_str());
}

void update_matrix(TRMTX T, double* X, int index)
{
    T[index][0] = X[0];
    T[index][1] = X[1];
    T[index][2] = X[2];
    T[3][index] = X[3];
}

void print_data(string fixed, string chf, string mobile, string chm, TRMTX T)
{
    cerr << "FIXED: " << fixed << endl;
    cerr << "MOBILE: " << mobile << endl;
    cerr << "CHAIN: " << chf << endl;
    cerr << "CHAIN: " << chm << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<3;j++)
        {
            cerr << T[i][j] << "   ";
        }
        cerr << endl;
    }
}

void read_cefile(string infile, string& fixed, string& chf, string& mobile, string& chm, TRMTX T)
{
    ifstream ifs;
    ifs.open(infile.c_str());
    if(!ifs)
    {
        cerr << "Could not open file: " << infile << endl;
        exit(EXIT_FAILURE);
    }

    int lc = 0;

    double X[4];
    string line;
    
    while(getline(ifs, line))
    {
        if(lc == 5)
        {
            get_filename_chain(line, fixed, chf);
        }
        if(lc == 6)
        {
            get_filename_chain(line, mobile, chm);
        }
        if(line.find("X2") != string::npos)
        {
            get_coords(trimmed(line), X);
            update_matrix(T, X, 0);
        }
        if(line.find("Y2") != string::npos)
        {
            get_coords(trimmed(line), X);
            update_matrix(T, X, 1);
        }
        if(line.find("Z2") != string::npos)
        {
            get_coords(trimmed(line), X);
            update_matrix(T, X, 2);
        }
        lc++;
    }
    ifs.close();

    //print_data(fixed, chf, mobile, chm, T);
}


