#include "read.h"

#include <fstream>
#include <iomanip>

//const string LOCN = "/var/www/3d-surfer/cgi-bin/tmp/";


void write_pdb(pdb& P, string ch)
{
    for(size_t i=0;i<P.atoms.size();i++)
    {
        cout << "ATOM  "
            << right
            << setw(5) << P.atoms[i].anum
            << "  "
            << left
            << setw(4) << P.atoms[i].atype
            << setw(3) << P.atoms[i].residue
            << " "
            //<< setw(1) << P.atoms[i].chain_id
            << setw(1) << ch
            //<< " "
            << right
            << setw(4) << P.atoms[i].rnum
            << "    "
            << fixed << setprecision(3)
            << setw(8) << P.atoms[i].axyz[0]
            << setw(8) << P.atoms[i].axyz[1]
            << setw(8) << P.atoms[i].axyz[2]
            << setw(6) << "1.00\n";
    }

}

void write_structures(pdb& pdb_F, pdb& pdb_R)
{
    //cerr << pdb_F.pname << " " << chf << " " << pdb_R.pname << " " << chm << endl;

    //string outfile = LOCN + pdb_F.pname + "_" + pdb_R.pname + "_rotated.pdb";
    //ofstream ofs;
    //ofs.open(outfile.c_str());
    //if(!ofs)
    //{
    //    cerr << "can't write pdb file: " << outfile << endl;
    //    exit(EXIT_FAILURE);
    //}

    write_pdb(pdb_F, "A");
    cout << "TER" << endl;
    write_pdb(pdb_R, "B");
    //ofs.close();
}

/**
 * point_transform:
 * @p: a Point.
 * @m: the Matrix representing the transformation to 
 * @origin: orgin of the reference frame
 * apply to the coordinates of @p.
 *
 * Transform the coordinates of @p according to @m. (p[] becomes m[][].p[]).
 */

void point_transform (double* p, TRMTX m)
{
    double x, y, z;
    x = m[0][0]*p[0] + m[0][1]*p[1] + m[0][2]*p[2] + m[3][0];
    y = m[1][0]*p[0] + m[1][1]*p[1] + m[1][2]*p[2] + m[3][1];
    z = m[2][0]*p[0] + m[2][1]*p[1] + m[2][2]*p[2] + m[3][2];
    p[0] = x; p[1] = y; p[2] = z;
}

// tranform coordinates for the target protein
void transform_coordinates(pdb& M, TRMTX T, pdb& R)
{
    vector<atom> V;
    for(size_t i=0;i<M.atoms.size();i++)
    {
        atom t = M.atoms[i];
        point_transform(t.axyz, T);
        V.push_back(t);
    }

    R.pname = M.pname;
    R.atoms = V;
}


void print_usage()
{
    cout << "\nUsage: Zorro ce_file \n" << endl;
    cout << "\nReport bugs to <vvenkatr@purdue.edu>\n" << endl;
}

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        print_usage();
        exit(EXIT_FAILURE);
    }
    if(argv[1] == NULL)
    {
        print_usage();
        exit(EXIT_FAILURE);
    }

    string infile(argv[1]);

    string fixed, mobile; // pdb files
    string chf, chm; // chains
    TRMTX T; // transformation matrix

    read_cefile(infile, fixed, chf, mobile, chm, T);

    //cerr << fixed << " " << mobile << endl;

    pdb pdb_F;
    read_protein(fixed, chf, pdb_F); // read first pdb
    //cerr << pdb_F.atoms.size() << endl;

    pdb pdb_M;
    read_protein(mobile, chm, pdb_M);// read second pdb
    //cerr << pdb_M.atoms.size() << endl;

    pdb pdb_R;
    transform_coordinates(pdb_M, T, pdb_R); // transform the atoms

    write_structures(pdb_F, pdb_R); // write the rotated structure

    // sayonara
    exit(EXIT_SUCCESS);
}
