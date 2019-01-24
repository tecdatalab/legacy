/*********************************************************************
*                           s i t u s . h                            *
**********************************************************************
* Header file for Situs C programs     URL: situs.biomachina.org     *
* (c) Willy Wriggers, Paul Boyle, and Pablo Chacon, 1998-2011        *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

/* widely used in most C programs */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <dirent.h>

#ifndef SITUS_H
#define SITUS_H

/* as of 2.7, used in colores, collage, eul2pdb, pdbsymm, lib_vio, and lib_eul: */
#define PI 3.14159265358979323846    

/* as of 2.7, used in colores, collage, eul2pdb, and lib_eul: */
#define ROT_CONV (PI/180.0) 

/* widely used in all PDB handling C programs */
typedef struct {  
	char  recd[7];        /*       1 -  6 */
	int   serial;         /*       7 - 11 */
	char  type[3];        /*      13 - 14 */
	char  loc[3];         /*      15 - 16 */
	char  alt[2];         /*           17 */
	char  res[5];         /*      18 - 21 */
	char  chain[2];       /*           22 */
	int   seq;            /*      23 - 26 */
	char  icode[2];       /*           27 */
	float x;              /*      31 - 38 */
	float y;              /*      39 - 46 */
	float z;              /*      47 - 54 */
	float occupancy;      /*      55 - 60 */
	float beta;           /*      61 - 66 */
	int   footnote;       /*      68 - 70 */
	char  segid[5];       /*      73 - 76 */
	char  element[3];     /*      77 - 78 */
	char  charge[3];      /*      79 - 80 */
	float weight;         /* mass of atom */
} PDB;  

#endif

/********************************************************************/
/*****  Old definitions for qpdb/qvol programs only          ********/
/********************************************************************/


#define SMAXS 100000  /* # of neural gas iteration steps            */ 
#define MAXCYCLE 8    /* max # of cycles in cluster analysis        */
#define NNMIN 2       /* minimum possible # of codebook vectors     */
#define NNMAX 50      /* maximum possible # of codebook vectors     */
#define MAXPDB 100000 /* maximum number of lines in pdb file        */
#define MAXDT (NNMAX*MAXCYCLE) 

typedef double Rseq3 [3];
typedef Rseq3 Rmat3NN [NNMAX];
typedef Rseq3 Rmat3DT [MAXDT];
typedef int IseqNN [NNMAX];
typedef double RseqNN [NNMAX];
typedef int ImatNNDT [NNMAX][MAXDT];

