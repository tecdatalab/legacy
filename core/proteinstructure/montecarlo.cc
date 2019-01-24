#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <cassert>

using namespace std;

#include <ANN/ANN.h>
#include "MersenneTwister.h"
#include "nr.h"
#include "charmm.h"
#include "atom.h"
#include "pdb.h"
#include "mdgraph.h"
#include "scoring.h"
#include "utils.h"
#include "montecarlo_options.h"

#define RANDOM_TRANSLATION_TRANSFORM 0.1 //0.5 for testing
#define RANDOM_ROTATION_TRANSFORM 0.05 //0.2 for testing

const double k_B = 1.38066e-23;
const double TEMPARATURE = 298.15; // 400 for testing
const double BETA = 1./(k_B * TEMPARATURE);
const double Z = 2.0 * M_PI * sqrt(M_PI);// constant used in Gaussian solvent exclusion calculation
const double DZ = 1./(M_PI*sqrt(M_PI)); // constant used in Gaussian solvent exclusion calculation for derivatives
const int NDIM = 6;
const DP GTOL=1.0e-06;
const double minimization_threshold = 15.0; /* score unit */
const char kBlankChars[] = " \t\n\r";

bool ORIGINAL = false;

const double WT_VDW = 0.2;
const double WT_GS = 0.3;
const double WT_ACP = 0.5;



typedef double Rotation[3][3];

MTRand objrnd;
vector<atom> REC, LIG, NLIG;
int nfunc, ndfunc;
double center_of_rotation[3];


void matrix_mul(Rotation a, Rotation b, Rotation c)
{
    for(int i = 0;i<3;i++)
    {
        for(int j = 0;j < 3; j++ )
        {
            c[i][j] = b[i][0] * a[0][j] + b[i][1] * a[1][j] + b[i][2] * a[2][j];
        }
    }
}

void xrotation(Rotation R, double alpha)
{
    double cosA = cos(alpha);
    double sinA = sin(alpha);

    R[0][0] = 1.0;
    R[0][1] = 0.0;
    R[0][2] = 0.0;

    R[1][0] = 0.0;
    R[1][1] = cosA;
    R[1][2] = -sinA;

    R[2][0] = 0.0;
    R[2][1] = sinA;
    R[2][2] = cosA;
}

void yrotation(Rotation R, double alpha)
{
    double cosA = cos(alpha);
    double sinA = sin(alpha);

    R[0][0] = cosA;
    R[0][1] = 0.0;
    R[0][2] = sinA;

    R[1][0] = 0.0;
    R[1][1] = 1.0;
    R[1][2] = 0.0;

    R[2][0] = -sinA;
    R[2][1] = 0.0;
    R[2][2] = cosA;
}

void zrotation(Rotation R, double alpha)
{
    double cosA = cos(alpha);
    double sinA = sin(alpha);

    R[0][0] = cosA;
    R[0][1] = -sinA;
    R[0][2] = 0.0;

    R[1][0] = sinA;
    R[1][1] = cosA;
    R[1][2] = 0.0;

    R[2][0] = 0.0;
    R[2][1] = 0.0;
    R[2][2] = 1.0;
}

// Gaussian distributed random translation
// rotation angle is also randomly distributed
void get_random_transform(double* TR)
{
    // obtain translation vector
    TR[0] = objrnd.randNorm(0., 1.) * RANDOM_TRANSLATION_TRANSFORM;
    TR[1] = objrnd.randNorm(0., 1.) * RANDOM_TRANSLATION_TRANSFORM;
    TR[2] = objrnd.randNorm(0., 1.) * RANDOM_TRANSLATION_TRANSFORM;

    // rotation
    
    TR[3] = objrnd.randNorm(0., 1.) * RANDOM_ROTATION_TRANSFORM;
    TR[4] = objrnd.randNorm(0., 1.) * RANDOM_ROTATION_TRANSFORM;
    TR[5] = objrnd.randNorm(0., 1.) * RANDOM_ROTATION_TRANSFORM;
}

// The acceptance rule for the MC simulation.
bool accept_move(double curr, double prev)
{
    double diff = curr - prev;
    if(diff < 0.)
        return true;
    else
    {
        double e = exp(-diff*BETA);
        double r = objrnd.rand();
        if(e > r)
            return true;
    }
    return false;
}

//  y --> R*(x-P) + (PT)
//  apply a rotation matrix R centered at P, and a translation T = PT-P
//  to a generic vector x(3), 
void apply_transformation(double* TR, vector<atom>& X, vector<atom>& Y)
{
    Rotation Rx, Ry, Rz, Rxy, ROT;
    xrotation(Rx, TR[3]);
    yrotation(Ry, TR[4]);
    zrotation(Rz, TR[5]);
    
    matrix_mul(Rx, Ry, Rxy); // Rxy = Ry * Rx
    matrix_mul(Rxy, Rz, ROT); // R = Rz * Ry * Rx

    double c[3];
    for(int i=0;i<3;i++)
    {
        c[i] = center_of_rotation[i] + TR[i];
    }

    double d[3], e[3];

    for(size_t i=0;i<X.size();i++)
    {
        atom a = X[i];
        d[0] = a.axyz[0] - c[0];
        d[1] = a.axyz[1] - c[1];
        d[2] = a.axyz[2] - c[2];

        e[0] = ROT[0][0] * d[0] + ROT[0][1] * d[1] + ROT[0][2] * d[2];
        e[1] = ROT[1][0] * d[0] + ROT[1][1] * d[1] + ROT[1][2] * d[2];
        e[2] = ROT[2][0] * d[0] + ROT[2][1] * d[1] + ROT[2][2] * d[2];

        a.axyz[0] = e[0] + c[0];
        a.axyz[1] = e[1] + c[1];
        a.axyz[2] = e[2] + c[2];

        Y.push_back(a);
    }
}

double vdw_score()
{
    vector<atom> L;
    if(ORIGINAL)
        L = LIG;
    else
        L = NLIG;

    double sigma, r, z, eps, tmp;
    double t;
    double A = (1./pow(0.6,12.)) - 2. * (1./pow(0.6, 6.));
    double B = -12.*(1./pow(0.6, 13.)) + 12. * (1./pow(0.6, 7.));
    double f_i = 1.0, f_j = 1.0;

    double vdw_rep = 0., vdw_attr = 0.;

    for(size_t i=0;i<REC.size();i++)
    {
        if(REC[i].INTF == 0)
            continue;

        if(REC[i].atype.at(0) == 'H')
            f_i = 0.4;
        else
            f_i = 1.0;

        for(size_t j=0;j<L.size();j++)
        {
            if(L[j].INTF == 0)
                continue;

            r = get_distance(REC[i].axyz, L[j].axyz);
            if(r > 8.)
                continue;

            if(L[j].atype.at(0) == 'H')
                f_j = 0.4;
            else
                f_j = 1.0;

            sigma = f_i * REC[i].chpar.RMIN + f_j * L[j].chpar.RMIN;
            z = sigma/r;
            eps = sqrt(REC[i].chpar.EMIN * L[j].chpar.EMIN);
            tmp = (eps*(pow(z, 12.) - 2.*pow(z, 6.)));
            t = 0.6 * sigma;


            if(r > t)
            {
                if(tmp < 0.)
                    vdw_attr += tmp;
                else
                   vdw_rep += tmp;
            }
            else
            {
                tmp = eps*(A + ((r-t)*(B/sigma)));
                vdw_rep += tmp;
            }
        }
    }


    return (vdw_rep + 5.*vdw_attr);
}

// Lazaridus-karplus
double gaussian_solvation_score()
{
    vector<atom> L;
    if(ORIGINAL)
        L = LIG;
    else
        L = NLIG;

    double r, f, n1, n2, t;
    double solv = 0.;

    for(size_t i=0;i<REC.size();i++)
    {
        if(REC[i].atype.at(0) == 'H') continue;
        if(REC[i].INTF == 0)
            continue;
        if(!REC[i].chpar.valset)
            continue;

        for(size_t j=0;j<L.size();j++)
        {
            if(L[j].atype.at(0) == 'H') continue;
            if(L[j].INTF == 0)
                continue;
            if(!L[j].chpar.valset)
                continue;

            r = get_squared_distance(REC[i].axyz, L[j].axyz);
            if(r > 64.0)
                continue;// from ZRANK

            n1 = 0.; n2 = 0.;

            n1 = (REC[i].chpar.DGFREE)/(r*r*REC[i].chpar.LAMBDA);
            t = (r - REC[i].chpar.RMIN)/REC[i].chpar.LAMBDA;
            f = exp(-(t*t));
            f *= L[j].chpar.AVOL;
            n1 *= f;

            n2 = (L[j].chpar.DGFREE)/(r*r*L[j].chpar.LAMBDA);
            t = (r - L[j].chpar.RMIN)/L[j].chpar.LAMBDA;
            f = exp(-t*t);
            f *= (REC[i].chpar.AVOL);
            n2 *= f;

            solv += (n1 + n2)/Z;
        }
    }
    return solv;
}

double contact_potential()
{
    vector<atom> L;
    if(ORIGINAL)
        L = LIG;
    else
        L = NLIG;

    int N = (int) L.size();
    if(N == 0)
    {
        cerr << "0 size vector. Exiting ..." << endl;
        exit(EXIT_FAILURE); 
    }
    
    ANNpointArray dataPts;
    dataPts = annAllocPts(N, 3);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<3;j++)
            dataPts[i][j] = L[i].axyz[j];
    }
    ANNkd_tree *T = new ANNkd_tree(dataPts, N, 3);
    ANNpoint f = annAllocPt(3);


    int l = 0, m = 0;
    double eps = 0.;
    double acp = 0.;

    double r2 = 36.;

    for(size_t i=0;i<REC.size();i++) // for the fixed protein
    {
        if(REC[i].acp_type == 0) continue;
        if(REC[i].INTF == 0)
            continue;
        f[0] = REC[i].axyz[0];
        f[1] = REC[i].axyz[1];
        f[2] = REC[i].axyz[2];

        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;

        int nnode = T->annkFRSearch(f, r2, 0, nnIdx, dists, eps); // search for n of them.
        if(nnode <= 0)
            continue;
        
        nnIdx = new ANNidx[nnode];
        dists = new ANNdist[nnode];
        T->annkFRSearch(f, r2, nnode, nnIdx, dists, eps); // search for n of them.
        l = REC[i].acp_type - 1;
        for(int p1=0; p1<nnode; p1++)
        {
            int j = nnIdx[p1];
            if(L[j].INTF == 0)
                continue;
            if(L[j].acp_type == 0) continue;
            m = L[j].acp_type - 1;
            acp += ZHANG_ACP[l][m];
        }
        delete [] nnIdx;
        delete [] dists;
    }
    acp *= (1./21.); // converting to kcal/mol

    delete T;
    annDeallocPts(dataPts);
    annClose();

    return acp;
}

// called by minimization recipes as the fuction to be minimized in docking.
// apply TR matrix to current and score the new structure, return the full
// score as the function value.
DP func(Vec_I_DP &x)
{
    nfunc++;
    return WT_VDW * vdw_score() + WT_ACP * contact_potential() + WT_GS * gaussian_solvation_score();
}



void dfunc(Vec_I_DP &x, Vec_O_DP &df)
{
    double TR[6];
    for(int i=0;i<6;i++)
    {
        TR[i] = x[i];
        df[i] = 0.0;
    }
    ndfunc++;

    vector<atom> Y;
    apply_transformation(TR, NLIG, Y);
    NLIG = Y;
    Y.clear();

    double dE_dr = 0.; // derivative of energy with respect to atom-atom distance
    double dr_dx[3]; // derivative of atom-atom distance along each coord. axis
    double q[3]; // position of j-atom relative to center of rotation
    double dsq; // squared distance between atoms
    double one_r; // 1/distance
    double dxyz[3];

    double sigma, r, z, eps, tmp, vdw_attr = 0., vdw_rep = 0.;
    double t, dt = 0.;
    double f_i = 1.0, f_j = 1.0;

    // VDW
    for(size_t i=0;i<REC.size();i++)
    {
        if(REC[i].INTF == 0)
            continue;
        for(size_t j=0;j<NLIG.size();j++)
        {
            if(NLIG[j].INTF == 0)
                continue;

            dxyz[0] = NLIG[j].axyz[0] - REC[i].axyz[0]; // calculate ij vector
            dxyz[1] = NLIG[j].axyz[1] - REC[i].axyz[1];
            dxyz[2] = NLIG[j].axyz[2] - REC[i].axyz[2];

            dsq = (dxyz[0] * dxyz[0]) + (dxyz[1] * dxyz[1]) + (dxyz[2] * dxyz[2]);
            // calculate distance squared

            if (dsq >= 64.0) continue;
            r = sqrt(dsq);

            sigma = f_i * REC[i].chpar.RMIN + f_j * NLIG[j].chpar.RMIN;
            z = sigma/r;
            eps = sqrt(REC[i].chpar.EMIN * NLIG[j].chpar.EMIN);
            tmp = (eps*(pow(z, 12.) - 2.*pow(z, 6.)));
            t = 0.6 * sigma;

            if(r > t)
            {
                if(tmp < 0.)
                    vdw_attr = tmp;
                else
                    vdw_rep = tmp;
            }
            else
            {
                vdw_rep = 1.;
            }

            dE_dr = 5.*vdw_attr + vdw_rep;
            dE_dr *= WT_VDW;

            one_r = 1.0/r;
            dr_dx[0] = dxyz[0]*one_r;
            dr_dx[1] = dxyz[1]*one_r;
            dr_dx[2] = dxyz[2]*one_r;

            // calculate the distance from the pivot point
            q[0] = NLIG[j].axyz[0] - center_of_rotation[0];
            q[1] = NLIG[j].axyz[1] - center_of_rotation[1];
            q[2] = NLIG[j].axyz[2] - center_of_rotation[2];

            df[0] += dE_dr * dr_dx[0]; // dE_dx
            df[1] += dE_dr * dr_dx[1]; // dE_dy
            df[2] += dE_dr * dr_dx[2]; // dE_dz
            df[3] += dE_dr * (dr_dx[2]*q[1] - dr_dx[1]*q[2]); // dE_dx-rotation
            df[4] += dE_dr * (dr_dx[0]*q[2] - dr_dx[2]*q[0]); // dE_dy-rotation
            df[5] += dE_dr * (dr_dx[1]*q[0] - dr_dx[0]*q[1]); // dE_dz-rotation
        }
    }

    double n1 = 0., n2 = 0., f = 0.;

    // Gaussian solvent exclusion
    for(size_t i=0;i<REC.size();i++)
    {
        if(REC[i].atype.at(0) == 'H') continue;
        if(REC[i].INTF == 0)
            continue;
        if(!REC[i].chpar.valset)
            continue;

        for(size_t j=0;j<NLIG.size();j++)
        {
            if(NLIG[j].atype.at(0) == 'H') continue;
            if(NLIG[j].INTF == 0)
                continue;
            if(!NLIG[j].chpar.valset)
                continue;

            dxyz[0] = NLIG[j].axyz[0] - REC[i].axyz[0]; // calculate ij vector
            dxyz[1] = NLIG[j].axyz[1] - REC[i].axyz[1];
            dxyz[2] = NLIG[j].axyz[2] - REC[i].axyz[2];

            dsq = (dxyz[0] * dxyz[0]) + (dxyz[1] * dxyz[1]) + (dxyz[2] * dxyz[2]);
            // calculate distance squared

            if (dsq >= 64.0) continue;
            r = sqrt(dsq);
            one_r = 1.0/r;

            n1 = 0.; n2 = 0.;

            n1 = (REC[i].chpar.DGFREE)/(r*r*REC[i].chpar.LAMBDA);
            t = (r - REC[i].chpar.RMIN)/REC[i].chpar.LAMBDA;
            f = exp(-(t*t));
            f *= NLIG[j].chpar.AVOL;
            dt = (t/REC[i].chpar.LAMBDA) + one_r;
            n1 *= (f*dt);

            n2 = (NLIG[j].chpar.DGFREE)/(r*r*NLIG[j].chpar.LAMBDA);
            t = (r - NLIG[j].chpar.RMIN)/NLIG[j].chpar.LAMBDA;
            f = exp(-t*t);
            f *= (REC[i].chpar.AVOL);
            dt = (t/REC[i].chpar.LAMBDA) + one_r;
            n2 *= (f*dt);

            dE_dr = (n1 + n2)/DZ;
            dE_dr *= WT_GS;
            
            dr_dx[0] = dxyz[0]*one_r;
            dr_dx[1] = dxyz[1]*one_r;
            dr_dx[2] = dxyz[2]*one_r;

            // calculate the distance from the pivot point
            q[0] = NLIG[j].axyz[0] - center_of_rotation[0];
            q[1] = NLIG[j].axyz[1] - center_of_rotation[1];
            q[2] = NLIG[j].axyz[2] - center_of_rotation[2];

            df[0] += dE_dr * dr_dx[0]; // dE_dx
            df[1] += dE_dr * dr_dx[1]; // dE_dy
            df[2] += dE_dr * dr_dx[2]; // dE_dz
            df[3] += dE_dr * (dr_dx[2]*q[1] - dr_dx[1]*q[2]); // dE_dx-rotation
            df[4] += dE_dr * (dr_dx[0]*q[2] - dr_dx[2]*q[0]); // dE_dy-rotation
            df[5] += dE_dr * (dr_dx[1]*q[0] - dr_dx[0]*q[1]); // dE_dz-rotation
        }
    }
}

void energy_minimize(double& energy_new, double* TR)
{
    int iter;
    DP fret;

    nfunc = ndfunc = 0;
    
    cerr << fixed << setprecision(4);
    //cerr << "Starting vector: (" << setw(8);
    //for(int i=0;i<6;i++)
    //    cerr << TR[i] << ",";
    //cerr << ")" << endl;
    //cerr << fixed << setprecision(6);

    Vec_DP p(NDIM);
    for(int i=0;i<6;i++)
        p[i] = TR[i];

    NR::dfpmin(p, GTOL, iter, fret, func, dfunc);
    
    //cerr << "Iterations: " << iter << endl;
    //cerr << "Func. evals: " << nfunc << endl;
    //cerr << "Deriv. evals: " << ndfunc << endl;

    //cerr << "Solution vector: (" << setw(8);
    //for(int i=0;i<6;i++)
    //    cerr << p[i] << ",";
    //cerr << ")" << endl;
    //cerr << "Func. value at solution " << setw(14) << fret << endl;
    energy_new = fret;
}


// the center of the interface is an appropriate point to rotate
// about when minimizing. this function finds that point, using CAs
// determine the interface residue list
// calculate the center of mass of interface residues (CA only)
void calculate_interface_centroid()
{
    int N = 0;
    for(int i=0;i<3;i++)
        center_of_rotation[i] = 0.;

    for(size_t i=0;i<LIG.size();i++)
    {
        if(LIG[i].INTF == 1)
        {
            if(LIG[i].atype.at(0) == 'H') continue;
            for(int j=0;j<3;j++)
                center_of_rotation[j] += LIG[i].axyz[j];
            N++;
        }
    }
    for(size_t i=0;i<REC.size();i++)
    {
        if(REC[i].INTF == 1)
        {
            if(REC[i].atype.at(0) == 'H') continue;
            for(int j=0;j<3;j++)
                center_of_rotation[j] += REC[i].axyz[j];
            N++;
        }
    }

    for(int j=0;j<3;j++)
        center_of_rotation[j] /= N;
}

int get_residue_CA(vector<atom>& S, int index)
{
    int f = 0;
    for(size_t i=0;i<S.size();i++)
    {
        if(S[i].residue == S[index].residue && S[i].atype == "CA" && S[i].rnum == S[index].rnum)
        {
            f = i;
            break;
        }
    }
    return f;
}

void set_interface_residues(vector<atom>& S, int rnum)
{
    for(size_t i=0;i<S.size();i++)
    {
        if(S[i].rnum == rnum)
            S[i].INTF = 1;
    }
}

void reset_residues()
{
    for(size_t i=0;i<REC.size();i++)
    {
        REC[i].INTF = -1;
    }
    for(size_t i=0;i<LIG.size();i++)
    {
        LIG[i].INTF = -1;
    }
}


void detect_interface_residues()
{
    int N = (int) REC.size();
    ANNpointArray dataPts;
    dataPts = annAllocPts(N, 3);


    for(int i=0; i<N; i++)
    {
        for(int j=0;j<3;j++)
            dataPts[i][j] = REC[i].axyz[j];
    }
    ANNkd_tree *TR = new ANNkd_tree(dataPts, N, 3);

    ANNpoint f = annAllocPt(3);

    for(size_t i=0;i<LIG.size();i++)
    {

        int ca = get_residue_CA(LIG, i);
        if(LIG[ca].INTF == 1)
            continue;

        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;

        f[0] = LIG[i].axyz[0];
        f[1] = LIG[i].axyz[1];
        f[2] = LIG[i].axyz[2];

        int nnode = TR->annkFRSearch(f, 100., 0, nnIdx, dists, 0.); // search for n of them.
        if(nnode <= 0)
            continue;
        
        set_interface_residues(LIG, LIG[i].rnum);
        delete [] nnIdx;
        delete [] dists;
    }

    delete TR;
    annDeallocPts(dataPts);
    

    N = (int) LIG.size();
    ANNpointArray BdataPts = annAllocPts(N, 3);
    for(int i=0; i<N; i++)
    {
        for(int j=0;j<3;j++)
            BdataPts[i][j] = LIG[i].axyz[j];
    }
    ANNkd_tree *TL = new ANNkd_tree(BdataPts, N, 3);

    
    
    for(size_t i=0;i<REC.size();i++)
    {
        int ca = get_residue_CA(REC, i);
        if(REC[ca].INTF == 1)
            continue;

        f[0] = REC[i].axyz[0];
        f[1] = REC[i].axyz[1];
        f[2] = REC[i].axyz[2];

        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;


        int nnode = TL->annkFRSearch(f, 100., 0, nnIdx, dists, 0.); // search for n of them.
        if(nnode <= 0)
            continue;

        set_interface_residues(REC, REC[i].rnum);

        delete [] nnIdx;
        delete [] dists;

    }

    delete TL;
    annDeallocPts(BdataPts);
    

    annClose();
}


void print_coords(string pdb_file_name)
{
    FILE *pdb_file ;
    if((pdb_file = fopen( pdb_file_name.c_str(), "w")) == NULL)
    {
        printf( "This file could not be opened.\nDying\n\n" ) ;
        exit(EXIT_FAILURE);
    }
    
    string s = REC[0].chain;
    
    for(size_t i=0;i<REC.size();i++)
    {
        if(s != REC[i].chain)
        {
            fputs("TER\n", pdb_file);
            s = REC[i].chain;
        }
        fprintf(pdb_file, "ATOM  %5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f\n", (int)i+1, REC[i].atype.c_str(), REC[i].residue.c_str(), REC[i].chain.c_str(), REC[i].rnum, REC[i].axyz[0], REC[i].axyz[1], REC[i].axyz[2]);        
    }
    fputs("TER\n", pdb_file);
    for(size_t i=0;i<LIG.size();i++)
    {
        fprintf(pdb_file, "ATOM  %5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f\n",
            (int)i+1, LIG[i].atype.c_str(), LIG[i].residue.c_str(), LIG[i].chain.c_str(), LIG[i].rnum, LIG[i].axyz[0], LIG[i].axyz[1], LIG[i].axyz[2]);
    }
    fclose(pdb_file);
}


// translate/rotate one partner in a docking complex
// this function minimizes the current position, and automatically
// accepts that position, since the score will be lower
// the relative orientation of two docking partners is represented by
// TR(6). Given a TR, a docking conformation can be defined and a score
// can be evaluated by the scorefxn. Here, we are minimizing the scorefxn
// by finding the best TR(6).


int main(int argc, char** argv)
{
	montecarlo_options options(argc,argv);
	if(!options.parse_successful())
	{
		cerr << "Usage: multilzerd_refine --input <ga file> --trials <trials_amount> --prefix <output_file_prefix> [-d <input_directory>] [-n <number_of_predictions]" << endl;
		exit(EXIT_FAILURE);
	}


	string ga_file = options.get_input_file();
	string directory = options.get_input_directory();
	string prefix = options.get_prefix();
	int results_number = options.get_results_number();
	int trials = options.get_trials();
	decoy_program_t decoy_program_type;

	// Read the information in the file
	ifstream ga_file_stream(ga_file.c_str(), ifstream::in);
	string main_complex_name, chains_arg, discard;
	int generations, population_size;
	double clash_threshold;

	string decoy_program;
	string score_type;
	// do it twice because we have "PDBID 1VCB", for example. The first one read PDBID to discard it
	ga_file_stream >> discard;
	ga_file_stream >> main_complex_name;
	// twice for the same reason
	ga_file_stream >> discard;
	ga_file_stream >> chains_arg;
	vector<string> chain_ids = split(chains_arg, ',');
	// Number of generations
	ga_file_stream >> discard;
	ga_file_stream >> generations;
	// Population size
	ga_file_stream >> discard;
	ga_file_stream >> population_size;
	// Clash threshold
	ga_file_stream >> discard;
	ga_file_stream >> clash_threshold;
	// Score type
	ga_file_stream >> discard;
	ga_file_stream >> score_type;
	// Decoy program
	ga_file_stream >> discard;
	ga_file_stream >> decoy_program;

	decoy_program_type = decoy_program.compare(DECOYPROGRAM_LZERD) == 0 ?
			decoy_program_lzerd : decoy_program_zdock;

	// load the pdb and prediction information (false, because we don't load hydrogens)
	vector<pdb> pdbs = load_pdbs(main_complex_name, chain_ids, directory, false);
	vector< vector<transformations*> > predictions = load_predictions(chain_ids, directory, decoy_program_type);
	vector<md_edge_set> population = load_population(ga_file_stream);

	// close the ga file because we don't need it anymore
	ga_file_stream.close();
	
	// This is used as a tmp holder to create the file for each optimized complex
	vector<vector<atom> > optimized;

	// Go through each of the predictions and run montecarlo optimization for each prediction
	size_t n = results_number == -1 ? population.size() : results_number;
	for(size_t population_index = 0; population_index < n; population_index++)
	{
		cout << "\nAnalyzing prediction #" << population_index + 1 << endl << endl;
		// Get the transformed atoms for the current prediction
		md_edge_set current_individual = population[population_index];
		vector< vector<atom> > transformed_proteins = current_individual.get_transformed_atoms(pdbs, predictions);

		// clear the previously used "receptor"/"ligand" vectors in order to add the new atoms
		REC.clear();
		LIG.clear();
		NLIG.clear();
		// Since the montecarlo optimization assumes that there are only 2 bodies (receptor/ligand)
		// we flip a coin for each chain in transformed_proteins to determine if it should be added
		// to the receptor or ligand (making sure that they both have at least one chain
		for(size_t chain_index = 0; chain_index < transformed_proteins.size(); chain_index++)
		{
			vector<atom>::iterator current_chain_begin = transformed_proteins[chain_index].begin();
			vector<atom>::iterator current_chain_end = transformed_proteins[chain_index].end();
			// test if it's the last one and one of the two bodies is empty (in which case we would just
			// add it to the empty one)
			if(chain_index == (transformed_proteins.size() - 1) && REC.empty())
			{
				REC.insert(REC.begin(), current_chain_begin, current_chain_end);

			}
			else if(chain_index == (transformed_proteins.size() - 1) && LIG.empty())
			{
				LIG.insert(LIG.begin(), current_chain_begin, current_chain_end);
			}
			else // normal case where we either have several chains to process or we are at the last one, but both REC
			// and LIG already have something in them. Just add it with uniform probability to either of them
			{
				if(objrnd.rand() > 0.5)
				{
					REC.insert(REC.begin(), current_chain_begin, current_chain_end);
				}
				else
				{
					LIG.insert(LIG.begin(), current_chain_begin, current_chain_end);
				}
			}
		}


		// get interface residues 10A
		detect_interface_residues();
		calculate_interface_centroid();

		double energy_new = 0.;
		ORIGINAL = true;
		double best_score = WT_VDW * vdw_score() + WT_ACP * contact_potential() + WT_GS * gaussian_solvation_score();

		ORIGINAL = false;

		double TR[6];
		for(int i=0;i<trials;i++)
		{
			// provide a random translation and rotation
			get_random_transform(TR);
			apply_transformation(TR, LIG, NLIG);

			energy_minimize(energy_new, TR);
			if(energy_new-best_score < minimization_threshold)
			{
				if(accept_move(energy_new, best_score))
				{
					cout << "Move accepted: " << best_score << " -> " << energy_new << endl;
					best_score = energy_new;
					LIG = NLIG;
					// reset receptor and ligand interface residues
					reset_residues();
					detect_interface_residues();
					calculate_interface_centroid();
				}
			}
			NLIG.clear();
// TODO: just the trace, remove it
/*
if(!(i % 10)) {
char numstring[10];
sprintf(numstring, "%05d", i + 1);
optimized.clear();
optimized.push_back(REC);
optimized.push_back(LIG);
write_complex(prefix + "t" + string(numstring), population_index, optimized);
}*/
		}
		cout << "Energy of Complex: " << best_score << endl;

		// create the file for the optimized complex
		optimized.clear();
		optimized.push_back(REC);
		optimized.push_back(LIG);
		write_complex(prefix, population_index, optimized);
	}

	exit(EXIT_SUCCESS);
}


