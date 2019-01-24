#include "contact.h"

int get_CA(string res, int rnum, string chain, vector<atom>& B)
{
    for(size_t i=0;i<B.size();i++)
    {
        if(B[i].atype == "CA")
        {
            if(res == B[i].residue && rnum == B[i].rnum && chain == B[i].chain)
                return (int) i;
        }
    }
    return -1;
}


void set_CA(map<int, int>& M, vector<atom>& A, vector<atom>& B, int& n)
{
    n = 0;
    int k = -1;
    for(size_t i=0;i<A.size();i++)
    {
        if(A[i].atype.compare("CA") == 0 && A[i].INTFA == 1)
        {
            k = get_CA(A[i].residue, A[i].rnum, A[i].chain, B);
            M.insert(make_pair((int) i, k));
            if(k != -1) n++;
        }
    }
}




// while looking for the c-alpha it will also set to 1 the INTF flag
// for all the residues that are on the interface.
// INTFA will only identify those residues on the interface and that are C-alpha at the same time
void set_interface_calpha(vector<atom>& X, string chain, string res, int rnum)
{
	for(size_t i=0;i<X.size();i++)
	{
		// all the atoms belonging to the residue defined by chain/res/rnum
		// will be marked as INTF=1
		if(res == X[i].residue && rnum == X[i].rnum && chain == X[i].chain)
		{
			X[i].INTF = 1;
			// and also, to identify an atom as a c-alpha interface one,
			// the INTFA flag is set to 1
			if(X[i].atype == "CA")
			{
				X[i].INTFA = 1;
			}
		}
	}
}

void get_calpha_pos(vector<atom>& X, ANNkd_tree* T)
{
    ANNpoint f = annAllocPt(3);
    double r2 = 100.0;
    for(size_t i=0;i<X.size();i++)
    {
        f[0] = X[i].axyz[0];
        f[1] = X[i].axyz[1];
        f[2] = X[i].axyz[2];
        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;

        int nnode = T->annkFRSearch(f, r2, 0, nnIdx, dists, 0.); 
        // set the CA atom of the residue to which this atom belongs            
        if(nnode > 0)
	{
            set_interface_calpha(X, X[i].chain, X[i].residue, X[i].rnum);
	}
    }
}

void check_CA(map<int, int>& M, int RL)
{
    int n1 = 0, n2 = 0;
    map<int, int>::iterator it;
    int f1, f2;
    for(it = M.begin();it!=M.end();++it)
    {
        f1 = it->first;
        f2 = it->second;
        if(f1 != -1)
            n1++;
        if(f2 != -1)
            n2++;
    }
    if(n2 == 0)
    {
        if(RL == 0)
        {
            cerr << "No matching CA in the 2 receptor structures" << endl;
        }
        else
        {
            cerr << "No matching CA in the 2 ligand structures" << endl;
        }
        cerr << "Original: " << n1 << endl;
        cerr << "Predicted: " << n2 << endl;
        exit(EXIT_FAILURE);
    }
    else if(n1 != n2)
    {
        cerr << "Some residues were not matched" << endl;
        cerr << "Original: " << n1 << endl;
        cerr << "Predicted: " << n2 << endl;
    }
}


void get_interface_residues(vector<atom>& R, vector<atom>& L)
{
    // get KDTree for R and L
    // get KDTree for receptor
    int N = (int) R.size();
    ANNpointArray RdataPts;
    RdataPts = annAllocPts(N, 3);

    for(int i=0; i<N; i++)
    {
        for(int j=0;j<3;j++)
            RdataPts[i][j] = R[i].axyz[j];
    }
    ANNkd_tree *TR = new ANNkd_tree(RdataPts, N, 3);

    N = (int) L.size();
    ANNpointArray LdataPts;
    LdataPts = annAllocPts(N, 3);

    for(int i=0; i<N; i++)
    {
        for(int j=0;j<3;j++)
	{
            LdataPts[i][j] = L[i].axyz[j];
	}
    }
    ANNkd_tree *TL = new ANNkd_tree(LdataPts, N, 3);

    get_calpha_pos(R, TL);
    get_calpha_pos(L, TR);

    // kdtree deletion
    delete TR;
    annDeallocPts(RdataPts);
    annClose();    
    // kdtree deletion
    delete TL;
    annDeallocPts(LdataPts);
    annClose();
}

// This function marks all atoms that belong to the interface between any pair of chains
// in a multi-chain complex
void get_interface_residues(vector< vector<atom> >& proteins)
{
	// If we have N proteins then the combinations are:
	// 1-2, 1-3, 1-4,....,1-N
	// then 2-3, 2-4, 2-5, ....,2-N
	// This double loop will go through all those combinations
	for(size_t receptor_index = 0; receptor_index < proteins.size(); receptor_index++)
	{
		for(size_t ligand_index = receptor_index + 1; ligand_index < proteins.size(); ligand_index++)
		{
			get_interface_residues(proteins[receptor_index], proteins[ligand_index]);
		}
	}
}



// LRMSD

void get_residue_list(vector<atom>& L,  vector<string>& IRES)
{
    set<string> VRES;
    for(size_t i=0;i<L.size();i++)
    {
        string str =  L[i].residue + int_to_string(L[i].rnum);
        VRES.insert(str);
    }
    set<string>::iterator it;
    for(it=VRES.begin();it!=VRES.end();it++)
    {
        IRES.push_back(*it);
    }
}

int get_CA(string res, vector<atom>& B)
{
    string str = "";
    for(size_t i=0;i<B.size();i++)
    {
        if(B[i].atype.compare("CA") == 0)
        {
            str = B[i].residue + int_to_string(B[i].rnum);
            if(str.compare(res) == 0)
                return (int) i;
        }
    }
    return -1;
}

void set_ligand_CA(map<int, int>& M, vector<string>& S, vector<atom>& A, vector<atom>& B, int& n)
{
    n = 0;
    string str;
    int k = -1;
    vector<string>::iterator it;
    for(size_t i=0;i<A.size();i++)
    {
        if(A[i].atype.compare("CA") == 0)
        {
            str = A[i].residue + int_to_string(A[i].rnum);
            it = find(S.begin(), S.end(), str);
            if(it != S.end())
            {
                k = get_CA(str, B);
                if(k != -1)
                {
                    n++;
                    M.insert(make_pair((int) i, k));
                }
            }
        }
    }
}

void check_ligand_CA(map<int, int>& M)
{
    int n1 = 0, n2 = 0;
    map<int, int>::iterator it;
    int f1, f2;
    for(it = M.begin();it!=M.end();++it)
    {
        f1 = it->first;
        f2 = it->second;
        if(f1 != -1)
            n1++;
        if(f2 != -1)
            n2++;
    }
    if(n2 == 0)
    {
        cerr << "No matching CA in the 2 ligand structures" << endl;
        cerr << "Original: " << n1 << endl;
        cerr << "Predicted: " << n2 << endl;
        exit(EXIT_FAILURE);
    }
}


// FNAT CALCULATION

// 2 residues are considered as being in contact if at least one atom of one residue is
// within 5Å of an atoms of the other.
void identify_contact_residues(vector<atom>& R, vector<atom>& L, ANNkd_tree* T, vector<string>& MP)
{
    int dim = 3;
    ANNpoint f = annAllocPt(dim);

    double r2 = 25.0;
    set<string> Z;
    
    for(size_t i=0;i<L.size();i++)
    {
        f[0] = L[i].axyz[0];
        f[1] = L[i].axyz[1];
        f[2] = L[i].axyz[2];

        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;

        int nnode = T->annkFRSearch(f, r2, 0, nnIdx, dists, 0.); // search for n of them.
        if(nnode <= 0)
            continue;

        nnIdx = new ANNidx[nnode];
        dists = new ANNdist[nnode];
        T->annkFRSearch(f, r2, nnode, nnIdx, dists, 0.); // search for n of them.

        for(int p1=0;p1<nnode;p1++)
        {
            int k = nnIdx[p1];
            string ar = R[k].chain + int_to_string(R[k].rnum) + R[k].residue;
            string br = L[i].chain + int_to_string(L[i].rnum) + L[i].residue;
            Z.insert(ar + "_" + br);
        }
        delete [] nnIdx;
        delete [] dists;

    }
    set<string>::iterator iter;
    for(iter=Z.begin();iter!=Z.end();++iter)
    {
        MP.push_back(*iter);
    }
}

// This is the multiple chain version of the previous function. Instead of assuming that we
// only have a receptor-ligand pair, we assume that we could have any number of structures
// For each pairwise combination of chains the pairwise identify_contact_residues will be
// called, adding the contacts found to the returned value MP.
// It will also return the contacts for each chain in contacts_per_chain; each position of this vector will represent one chain
void identify_contact_residues(vector< vector<atom> > proteins, vector<string>& MP, vector< vector<string> >& contacts_per_chain)
{
	// reserve the space for the contacts info, for each chain
	contacts_per_chain.resize(proteins.size());
	// If we have N proteins then the combinations are:
	// 1-2, 1-3, 1-4,....,1-N
	// then 2-3, 2-4, 2-5, ....,2-N
	// This double loop will go through all those combinations
	for(size_t receptor_index = 0; receptor_index < proteins.size(); receptor_index++)
	{
		for(size_t ligand_index = receptor_index + 1; ligand_index < proteins.size(); ligand_index++)
		{
			vector<string> current_chain_contacts;
			vector<atom> receptor = proteins[receptor_index];
			vector<atom> ligand = proteins[ligand_index];
			int receptor_size = (int) receptor.size();
			ANNpointArray AdataPts;
			AdataPts = annAllocPts(receptor_size, 3);

			for(int i=0; i<receptor_size; i++)
			{
				for(int j=0;j<3;j++)
				{
					AdataPts[i][j] = receptor[i].axyz[j];
				}
			}
			ANNkd_tree *TR = new ANNkd_tree(AdataPts, receptor_size, 3);

			identify_contact_residues(receptor, ligand, TR, current_chain_contacts);

			MP.insert(MP.end(), current_chain_contacts.begin(), current_chain_contacts.end());
			contacts_per_chain[receptor_index].insert(contacts_per_chain[receptor_index].end(),
								current_chain_contacts.begin(), current_chain_contacts.end());
			contacts_per_chain[ligand_index].insert(contacts_per_chain[ligand_index].end(),
								current_chain_contacts.begin(), current_chain_contacts.end());

			// kdtree deletion
			delete TR;
			annDeallocPts(AdataPts);
		}
	}
}

// Method used to help on the calculation of RMSD. It returns a MAP where the key is a pair of predicted_chain/original_chain (the indices)
// and the value of the map is a vector with the number of all the atoms that match between the original chain and the predicted one
int get_matching_atoms(vector<vector<atom> >& A, vector<vector<atom> >& B, map<pair<int,int>, vector<int> >& M, bool interface_only)
{
	if(A.size() != B.size())
	{
		cerr << "Number of chains in protein are not the same" << endl;
		exit(EXIT_FAILURE);
	}

	int n = 0, d = -1;
	for(size_t i=0;i<A.size();i++)
	{
		vector<atom> F = A[i];
//cout << "A " << F[i].rnum << " " << F[i].chain << " " << F[i].residue << "\n";
		vector<int> V;
		d = -1;
		for(size_t j=0;j<B.size();j++)
		{
			vector<atom> G = B[j];
//cout << "B " << G[j].rnum << " " << G[j].chain << " " << G[j].residue << "\n";
			if(F[0].chain != G[0].chain)
			{
				continue;
			}
			else
			{
				d = (int) j;
				for(size_t k=0;k<F.size();k++)
				{
/*cout << "C " << F[k].rnum << " " << F[k].residue << " " << F[k].atype << " " <<
	G[k].rnum << " " << G[k].residue << " " << G[k].atype <<"\n";
cout << (F[k].residue.compare(G[k].residue)) << " " << (F[k].rnum != G[k].rnum) << " " <<
	(F[k].atype.compare("CA") != 0) << " " << (G[k].atype.compare("CA") != 0) << endl;
*/					// prune all those residues that don't match and only keep C-alpha
					if((F[k].residue.compare(G[k].residue)) || (F[k].rnum != G[k].rnum)
						|| (F[k].atype.compare("CA") != 0) || (G[k].atype.compare("CA") != 0))
					{
//cout << "D breaking\n";
						continue;
					}
					// if interface_only then only process those
					// only one will be marked as interface
					if(interface_only && (F[k].INTFA == -1) && (G[k].INTFA == -1))
					{
//cout << "E breaking\n";
						continue;
					}
//cout << "F added\n";
					V.push_back(k);
				}
			}
		}
		//if(V.size() == 0)
		//{
			//cerr << "Error: could not find matching CA" << endl;
			//exit(EXIT_FAILURE);
		//}
		if(V.size())
		{
			M.insert(make_pair(make_pair(i,d), V));
			n += (int) V.size();
		}
	}
	return n;
}
