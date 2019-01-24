#include "read.h"
#include <unistd.h>

bool is_member(vector<int>& X, int k)
{
    vector<int>::iterator pos;
    pos = find(X.begin(), X.end(), k);
    if(pos == X.end())
        return false;
    else
        return true;
}

void print_coordinates(pdb& P, vector<vector<int> >& V)
{
    for(size_t i=0;i<V.size();i++)
    {
        vector<int> X = V[i];
        string apName = P.pname;
        stringstream ss;
        ss << (int) (i+1);
        ss >> apName;
        string RES = "results_" + apName + ".txt";
        ofstream ofs(RES.c_str());
   
        vector<int> M;

        for(size_t j=0;j<P.atoms.size();j++)
        {
            if(is_member(X, P.atoms[j].rnum))
            {
                M.push_back(j);
            }
        }

        ofs << 3 << endl;
        ofs << M.size() << endl;

        for(size_t j=0;j<M.size();j++)
        {

            ofs << P.atoms[M[j]].axyz[0] << " "
                << P.atoms[M[j]].axyz[1] << " "
                << P.atoms[M[j]].axyz[2]
                << endl;
        }

        ofs.close();
    }

}

int main(int argc, char** argv)
{
    int c;
    extern char *optarg;
    extern int optind;
    string pdbfile = "", visfile = "";
    bool errflg = false;
    int N = 0;

    while((c = getopt(argc, argv, "hp:v:")) != EOF)
    {
        switch(c)
        {
        case 'p':
            pdbfile.assign(optarg);
            break;
        case 'v':
            visfile.assign(optarg);
            break;
        case 'h':
            errflg = true;
            break;
        case '?':
            errflg = true;
            break;
        }
    }

    if(pdbfile.length() == 0 || visfile.length() == 0)
    {
        cerr << "No file for processing" << endl;
        errflg = true;
    }

    if(errflg)
    {
        cerr << "Usage: getpoints [-p<pdbfile>] [-v<visgrid>] \n";
        cerr << "-p pdbfile      pdb file to be processed" << endl;
        cerr << "-v visgridfile  file containing Visgrid output" << endl;
        exit(EXIT_SUCCESS);
    }

    pdb P;
    read_protein(pdbfile, P);

    vector<vector<int> > res_list;
    read_visgrid_output(visfile, res_list);

    print_coordinates(P, res_list);

    exit(EXIT_SUCCESS);
}
