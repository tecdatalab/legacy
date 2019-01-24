/*
 * can be used interchangeably with ".pqr" files
 */

#include <string>
#include <vector>
using namespace std;


#ifndef _ATOM_H_
#define _ATOM_H_

class atom
{
    public:
        int         anum;   // id of the atom
        string      atype;  // type of the atom
        int         rnum;   // residue id
        double      axyz[3]; // xyz coordinates
        string      chain_id; // chain id
        string      residue; // name of residue
        string      atom_char;

        atom(int m_num, string m_type, int m_rid, double* m_xyz, string m_residue,
            string achar, string chain)
        {
            anum = m_num;
            rnum = m_rid;
            atype = m_type;
            for(int i=0;i<3;i++)
            {
                axyz[i] = m_xyz[i];
            }
            residue = m_residue;
            atom_char = achar;
            chain_id = chain;
        }

        atom(){}
        ~atom(){}
};

#endif

#ifndef _BOND_H_
#define _BOND_H_

class bond
{
    public:
        int bstart;
        int bend;
        bond(int mstart, int mend)
        {
            bstart = mstart;
            bend = mend;
        }
        bond(){}
        ~bond(){}
};

#endif


/**
 * PDB header file
 * stores atom coordinates,types and bond information
 *
 * @author vishwesh venkatraman
 * @date   19/08/2007
 */

#ifndef _PDB_H_
#define _PDB_H_

class pdb
{
public:
    string pname;
    vector<atom> atoms;
    vector<bond> bonds;
    pdb(string m_name, vector<atom>& m_atoms, vector<bond>& m_bonds)
    {
        pname = m_name;
        atoms = m_atoms;
        bonds = m_bonds;
    }
    pdb(string m_name, vector<atom>& m_atoms)
    {
        pname = m_name;
        atoms = m_atoms;
    }
    pdb(){}
    ~pdb(){};
};

#endif
