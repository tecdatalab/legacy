#ifndef _MARCHING_CUBES_SURFACE_H_
#define _MARCHING_CUBES_SURFACE_H_

#include "point_lattice.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>		   // Open Graphics Library (OpenGL) header
#include <GL/glut.h>	   // The GL Utility Toolkit (GLUT) Header
#endif

// This was originally a struct but it was turned into
// a class to have a default constructor that initializes
// all values to zero
class GLvector
{
  public:
	GLfloat fX;
	GLfloat fY;
	GLfloat fZ;
	GLvector() : fX(0), fY(0), fZ(0) {}
};

/*
 * Given a point_lattice instance and an isovalue,
 * this class can issue appropriate OpenGL commands to draw
 * a surface using marching cubes.
 *
 * The implementation is based on a Marching Cubes
 * Example Program by Cory Bloyd (corysama@yahoo.com).
 * The code is provided with a public domain intention.
 * For a description of the algorithm go to
 * http://astronomy.swin.edu.au/pbourke/modelling/polygonise/
 * or
 * http://undergraduate.csse.uwa.edu.au/units/CITS4241/Project/references/Lorensen-Cline-brief.pdf
 */
class marching_cubes_surface
{
  private:
	point_lattice* lattice;
	// If this is set to true, instead of drawing filled triangles
	// the draw method will only draw triangle edges
	bool wireframe_only;
	// a2fVertexOffset lists the positions, relative to vertex0,
	// of each of the 8 vertices of a cube
	static unsigned long a2fVertexOffset[8][3];
	// a2iEdgeConnection lists the index of the endpoint vertices
	// for each of the 12 edges of the cube
	static GLint a2iEdgeConnection[12][2];
	// a2fEdgeDirection lists the direction vector (vertex1-vertex0)
	// for each edge in the cube
	static GLfloat a2fEdgeDirection[12][3];
	// For any edge, if one vertex is inside of the surface and the other
	// is outside of the surface then the edge intersects the surface
	// For each of the 8 vertices of the cube can be two possible states :
	// 		either inside or outside of the surface
	// For any cube the are 2^8=256 possible sets of vertex states
	// This table lists the edges intersected by the surface for all
	// 256 possible vertex states
	// There are 12 edges.  For each entry in the table, if edge #n
	// is intersected, then bit #n is set to 1
	static GLint aiCubeEdgeFlags[256];
	//  For each of the possible vertex states listed in aiCubeEdgeFlags
	// there is a specific triangulation of the edge intersection points.
	// a2iTriangleConnectionTable lists all of them in the form of
	//  0-5 edge triples with the list terminated by the invalid value -1.
	//  For example: a2iTriangleConnectionTable[3] list the 2 triangles
	// formed when corner[0] and corner[1] are inside of the surface,
	// but the rest of the cube is not.
	static GLint a2iTriangleConnectionTable[256][16];
  public:
	// Just provide the lattice that will be drawn
	// Optionally we can indicate if we just want to draw
	// triangle edges (default is to draw filled triangles
	marching_cubes_surface(point_lattice* in_lattice,
							bool in_wireframe_only=false);
	inline void set_wireframe(bool in_wireframe_only) {
		wireframe_only = in_wireframe_only;
	}
	inline bool is_wireframe_only() {
		return wireframe_only;
	}
	// This method will be invoked by the main OpenGL coordinator
	// program. It just pushes the current stack, calls the
	// appropriate triangle drawing methods and then returns
	// The isovalue provided is used to determine what is in or
	// or out of the surface
	void draw(float isovalue);
	// Applies the marching cubes algorithm to the cube whose bottom-left-front
	// corner is located at fX, fY, fZ
	void march_cube(unsigned long fX, unsigned long fY, unsigned long fZ,
					float isovalue);
	// Finds the approximate point of intersection of the surface
	// between two points with the values fValue1 and fValue2
	GLfloat get_offset(GLfloat fValue1, GLfloat fValue2, GLfloat fValueDesired);
	// get_normal() finds the gradient of the scalar field at a point
	// This gradient can be used as a very accurate vertx normal for
	// lighting calculations
	void get_normal(GLvector &rfNormal, GLfloat fX, GLfloat fY, GLfloat fZ);
	GLvoid normalize_vector(GLvector &rfVectorResult, GLvector &rfVectorSource);
};

#endif
