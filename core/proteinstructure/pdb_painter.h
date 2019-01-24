#ifndef _PDB_PAINTER_H_
#define _PDB_PAINTER_H_

#include "pdb.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>		   // Open Graphics Library (OpenGL) header
#include <GL/glut.h>	   // The GL Utility Toolkit (GLUT) Header
#endif

/*
 * Draws a representation of a protein extracted from a pdb
 * instance. Initially only a c-alpha representation is intended.
 */
class pdb_painter
{
 private:
	pdb* protein;
 public:
  pdb_painter(pdb* in_protein);
  void draw_calpha();
};

#endif
