#ifndef __BOX_H__
#define __BOX_H__

#include <stdlib.h>
#include <string.h>
#include <GL/glut.h>
#include <malloc.h>
#include <memory.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <valarray>
#include <list>
#include <vector>
#include <fstream>
#include <iomanip>
#include "Grider.h"

#include "glui.h"
#define RANDOM ((float)rand())/RAND_MAX
//#define RETURN_FALSE m_boxAry.free(); m_polygons.free(); return false

using namespace std;

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#define GET_LINE f.getline(buffer, 200); while (strstr(buffer, "#")!=NULL) f.getline(buffer, 200);
#define pi 3.14159265359

unsigned int g_polygonNum;
extern int g_radiogroup, g_materials, g_material_faces, g_Pausing, g_viewOctree, g_FrustrumCulling, g_viewPDB, g_viewGrid, g_threshold;
extern float g_unit;
extern float g_minOctreeBoxSize;
extern float g_materialv[];

class CVector {
public:
	inline  CVector() { X = Y = Z = 0;}
	inline  CVector operator+(CVector Vecf) const { return CVector(X + Vecf.X, Y + Vecf.Y,Z+Vecf.Z); }
	inline  CVector operator-(CVector Vecf) const { return CVector(X - Vecf.X, Y - Vecf.Y,Z-Vecf.Z); }
	inline  CVector operator*(float f) const { return CVector(X*f, Y*f, Z*f); }
	inline  CVector(float initX, float initY, float initZ) { X=initX; Y=initY; Z=initZ;}
	inline  void operator/=(float f) { X/=f; Y/=f; Z/=f; }
	inline  float len() { return sqrt(X*X+Y*Y+Z*Z); }
	inline	CVector Max(CVector Vecf) { return CVector(max(X, Vecf.X), max(Y, Vecf.Y), max(Z, Vecf.Z)); }
	inline	CVector Min(CVector Vecf) { return CVector(min(X, Vecf.X), min(Y, Vecf.Y), min(Z, Vecf.Z)); }
	float X,Y,Z;
};

class CVertex {
public:
	CVertex() { hasNormal = 0; hasColor = 0; }
	CVector p, normal, color;
	unsigned char hasNormal:1;
	unsigned char hasColor:1;
};

class CBox {
public:
	inline void init(CVector pixel) {
		minV = maxV = pixel;
	}
	inline void mergePixel(int index, CVector pixel) { 
		this->index = index; 
		maxV = maxV.Max(pixel); minV = minV.Min(pixel);
	}
	int index; CVector minV, maxV;
};

class CPolygon : public valarray<CVertex> {
public:
	bool visible;

	inline CBox getBox(int index) {
		unsigned int i; CBox box;
		if (size()==0) return box;
		box.init((*this)[0].p);
		for (i=1; i<size(); i++) {
			box.mergePixel(index, (*this)[i].p);
		}
		return box;
	}
};

#define OFFSET(x, y, z) (x+y*m_nSize+z*m_nSize*m_nSize)

struct CubePixel { unsigned int count; bool bVisible; };

class OctTree {
public:
	OctTree() { for (int i=0; i<8; i++) m_child[i] = NULL; }
	OctTree* m_child[8];
	list<unsigned int> m_boxIndices;
	CVector m_minCorner, m_maxCorner;
};

#define GL_VERTEX_VEC(p) glVertex3f(p.X, p.Y, p.Z);

class CWorldCube {
public:
	void freeAll() {
	}
	void BuildWorld(valarray<ATOM>& atoms) {
		unsigned int i;
		CVector v;
		v.X = atoms[1].v.x; v.Y = atoms[1].v.y; v.Z = atoms[1].v.z;
		m_minCorner = m_maxCorner = v; 
		for (i=1; i<atoms.size(); i++) {	
			v.X = atoms[i].v.x; v.Y = atoms[i].v.y; v.Z = atoms[i].v.z;
			m_maxCorner = m_maxCorner.Max(v); m_minCorner = m_minCorner.Min(v); 
		}
		m_worldCenter = (m_maxCorner+m_minCorner); m_worldCenter /= 2.0f;
		m_fWorldSize = max(m_maxCorner.X-m_minCorner.X, m_maxCorner.Y-m_minCorner.Y); 
		m_fWorldSize = max(m_fWorldSize, m_maxCorner.Z-m_minCorner.Z);
		m_minCorner = m_maxCorner = m_worldCenter;
		float radius = m_fWorldSize/2;
		m_minCorner.X -= radius; m_minCorner.Y -= radius; m_minCorner.Z -= radius;
		m_maxCorner.X += radius; m_maxCorner.Y += radius; m_maxCorner.Z += radius;
	}

	static inline void drawBox(CVector minV, CVector maxV) {
		CVector corners[8];
		corners[0] = corners[1] = corners[2] = corners[3] = minV;
		corners[1].X = maxV.X; corners[2].Y = maxV.Y; corners[3].Z = maxV.Z;
		corners[4] = corners[5] = corners[6] = corners[7] = maxV;
		corners[4].X = minV.X; corners[5].Y = minV.Y; corners[6].Z = minV.Z;

#ifdef _USE_OPENGL_
		glBegin (GL_LINES);	glColor3f(1.0f, 1.0f, 1.0f);
		GL_VERTEX_VEC(corners[0]); GL_VERTEX_VEC(corners[1]);
		GL_VERTEX_VEC(corners[0]); GL_VERTEX_VEC(corners[2]);
		GL_VERTEX_VEC(corners[0]); GL_VERTEX_VEC(corners[3]);
		GL_VERTEX_VEC(corners[7]); GL_VERTEX_VEC(corners[4]);
		GL_VERTEX_VEC(corners[7]); GL_VERTEX_VEC(corners[5]);
		GL_VERTEX_VEC(corners[7]); GL_VERTEX_VEC(corners[6]);
		GL_VERTEX_VEC(corners[1]); GL_VERTEX_VEC(corners[5]);
		GL_VERTEX_VEC(corners[1]); GL_VERTEX_VEC(corners[6]);
		GL_VERTEX_VEC(corners[2]); GL_VERTEX_VEC(corners[4]);
		GL_VERTEX_VEC(corners[2]); GL_VERTEX_VEC(corners[6]);
		GL_VERTEX_VEC(corners[3]); GL_VERTEX_VEC(corners[4]);
		GL_VERTEX_VEC(corners[3]); GL_VERTEX_VEC(corners[5]);
		glEnd ();
#endif
	}

	//GLint m_viewport[4]; GLdouble m_mvmatrix[16], m_projmatrix[16];				// Transform Matrix
	float m_fStep, m_fWorldSize;
	float m_minSize;
	CVector m_minCorner, m_maxCorner, m_worldCenter;
};

class CWorld
{
public:
	CWorld() { init(); };
	void init() {
		cmdIndex = m_x0 = m_y0 = m_z0 = m_x1 = m_y1 = m_z1 = m_scale = m_ox = m_oy = m_oz = 0.0f;
		m_zNear = 0.01f; m_zFar = 999.0f;
		m_rx = m_ry = m_rz = m_f = m_vx = m_vy = m_vz = m_dx = m_dy = m_dz = m_ux = m_uy = m_uz = 0.0f;
		m_scale = 1.0f; m_f=60.0f; m_vz=-1.0f; 
		m_lx = 0.0f; m_ly = 2.0f; m_lz = 3.0f;
		m_ambientR = 0.7; m_ambientG = 0.5; m_ambientB = 0.3;
		m_diffuseR = 0.4; m_diffuseG = 0.2; m_diffuseB = 0.5;
		srand( (unsigned)time(NULL));
		g_minOctreeBoxSize = 10.0;
		g_viewOctree = 0;
		g_FrustrumCulling = 0;
		g_viewPDB = 1;
		g_viewGrid = 0;
		g_unit = 0.5;
	}

	~CWorld(void){  }

	void freeAll() { m_worldBox.freeAll();}

	bool loadfile(char* fname) {
		int i, j, k; string base;
		Vecf maxVec, minVec, vecCenter;
		int atomnum, gsize; 
		float maxspan, gspan, unit;

		CReader reader;
		if (!reader.read(fname) || reader.size()==0) {
			printf("Can't open %s \n", fname); return -2; 
		}

		for (atomnum=0, i=0; i<reader.size(); i++) {
			if (strncmp("ATOM", reader[i], 4) != 0 && strncmp("HETATM", reader[i], 6) != 0) continue;
			atomnum++;
		}

		m_atoms.resize(atomnum+1);
		for (i=0, k=1; i<reader.size(); i++) {
			if (strncmp("ATOM", reader[i], 4) != 0 && strncmp("HETATM", reader[i], 6) != 0) continue;
			_ATOM* p = (_ATOM*)reader[i];	ATOM& atom = m_atoms[k++];

			atom.bHetAtom = false;
			if (strncmp("HETATM", p->Atom, 6) == 0) atom.bHetAtom = true;
			sscanf(p->serial, "%d", &atom.serial);
			for (base = "", j=0; j<4; j++) {
				if (p->name[j] != ' ') base += p->name[j];
			}	strcpy(atom.name, CSTR(base));
			for (base = "", j=0; j<3; j++) {
				if (p->resName[j] != ' ') base += p->resName[j];
			}	strcpy(atom.resName, CSTR(base));
			atom.chainID = p->chainID;
			sscanf(p->resSeq, "%d", &atom.resSeq);
			atom.iCode = p->iCode;
			sscanf(p->x, "%f", &atom.v.x);
			sscanf(p->y, "%f", &atom.v.y);
			sscanf(p->z, "%f", &atom.v.z);
			atom.radius = get_radius(atom.resName, atom.name); 
			atom.binded = 0;
			atom.lv = 0;

			if (k==2) maxVec = minVec = atom.v;
			maxVec.x = max(atom.v.x, maxVec.x); maxVec.y = max(atom.v.y, maxVec.y); maxVec.z = max(atom.v.z, maxVec.z);
			minVec.x = min(atom.v.x, minVec.x); minVec.y = min(atom.v.y, minVec.y); minVec.z = min(atom.v.z, minVec.z);
		}

		vecCenter = Vecf((maxVec.x+minVec.x)/2,(maxVec.y+minVec.y)/2,(maxVec.z+minVec.z)/2);
		maxspan = (maxVec-minVec).span()/2;

		gspan = (float)(4 + 2*(maxspan+3.0));
		gsize = ((int)(gspan/g_unit+1)/2)*2;
		unit = gspan/gsize;
		Vecf gcenter(gsize/2, gsize/2, gsize/2);

		for (i=1; i<(int)m_atoms.size(); i++) {
			ATOM& atom = m_atoms[i];
			atom.gv = (atom.v - vecCenter)/(float)g_unit+gcenter;
		} 

		m_grid.resize(gsize, gsize, gsize);
		VERIFY(m_atoms.size() < 65535, "atoms counts exceed 65535!");

		m_fillCnt = 0;
		for (i=1; i<m_atoms.size(); i++) {
			ATOM& atom = m_atoms[i]; VOXEL& vox = m_grid.vox(atom.gv);
			if (atom.bHetAtom) continue;
			VERIFY(atom.radius>0, "atom.radius>0");
			float radius = 1.4+atom.radius; radius = radius/unit;
			vox.m_type = 'F'; vox.m_atomserial = i;

			CRange r(m_grid, atom.gv, radius);
			while (r.next()) { 
				VOXEL& vox = r.vox();
				float dist = (r.v-atom.gv).len();
				if (dist <= radius) { 
					if (vox.m_type != 'F') m_fillCnt++;
					vox.m_type = 'F';
					if (vox.m_atomserial != -1) {
						ATOM& atomX = m_atoms[vox.m_atomserial];
						float distX = (r.v-atomX.gv).len();
						if (dist/atom.radius > distX/atomX.radius)
							continue;
					} vox.m_atomserial = i;
				}
			}
		}

		int n, dim=gsize;
		FILE* fp = fopen("data.bin", "wb");
		float one = 1.0f, zero = 0.0f;
		n=0;
		for (i=0; i<dim; i++) {
			for (j=0; j<dim; j++) { 
				for (k=0; k<dim; k++) {
					if (m_grid[n++].m_type == 'F')
						fwrite(&one, sizeof(float), 1, fp); 
					else
						fwrite(&zero, sizeof(float), 1, fp); 
				}
			}
		}
		fclose(fp);

		if (m_seedNum < 0) {
			return true;
		}

		ATOM& fe = m_atoms[m_seedNum];

		int size = m_atoms.size()-1; CVector Vecf; init(); 

		m_worldBox.freeAll();
		m_worldBox.m_minSize = g_minOctreeBoxSize;
		m_worldBox.BuildWorld(m_atoms);

		center = m_worldBox.m_maxCorner + m_worldBox.m_minCorner; center/=2.0f;
		maxV = m_worldBox.m_maxCorner; minV = m_worldBox.m_minCorner;

		float objsize = max((maxV-center).len(), (center-minV).len());
		m_x0 = minV.X-objsize; m_y0 = minV.Y-objsize; m_z0 = minV.Z-objsize;
		m_x1 = maxV.X+objsize; m_y1 = maxV.Y+objsize; m_z1 = maxV.Z+objsize;
		m_scale = 1.0f;
		m_ox = -center.X; m_oy = -center.Y; m_oz = -center.Z;
		m_rx = 0.0f; m_ry = 0.0f; m_rz = 0.0f;
		m_f  = 60.0f;
		m_vx = 0.0f; m_vy = 0.0f, m_vz = 4*objsize;
		m_dx = 0.0f; m_dy = 0.0f; m_dz = -1.0f;
		m_ux = 0.0f; m_uy = 1.0f; m_uz = 0.0f;

		m_objSize=maxV-minV; m_objSize/=2.0f;
		m_ox += RANDOM*m_objSize.X; m_oy = RANDOM*m_objSize.Y; m_oz = RANDOM*m_objSize.Z;
		m_lx = m_ox; m_ly = m_oy; m_lz = m_oz;
		m_bPaused = true;

		calculate_matrix();
		return true;
	}

	int maxCnt(int radius) {
		int i, j, k, total=0;
		int r2 = radius*radius;
		for (i=-radius; i<=radius; i++) {
			for (j=-radius; j<=radius; j++) {
				for (k=-radius; k<=radius; k++) {
					if (i*i+j*j+k*k <= r2)
						total++;
				}
			}
		}
		return total-1;
	}

	void report() {
		ATOM& fe = m_atoms[m_seedNum];
		FILE *fp = fopen("out.txt","w");
		fprintf(fp, "index\tsize\tfilled\tmaxFilled\t1-filled/maxFilled\n");

		int i;
		for (i=1; i<m_grid.width; i++) {
			if (flood(i, fp))
				break;
			fflush(fp);
		} fclose(fp);
	}

	bool flood(int threshold, FILE* fp=NULL) {
		int i;
		ATOM& fe = m_atoms[m_seedNum];

		Vecf seed; bool seeded = false;
		for (i=1; i<m_grid.depth; i++) {
			CRange r(m_grid, fe.gv, i);
			while (r.next()) { 
				VOXEL& vox = r.vox();
				float dist = (r.v-fe.gv).len();
				if (dist <= i) {
					if (vox.m_type != 'F') { seed = r.v; seeded = true; break;}
				}
			}
			if (seeded) break;
		}

		int minCnt;
		minCnt = 0;

		int size = m_grid.depth*m_grid.width*m_grid.height;
		for (i=0; i<size; i++) {
			VOXEL& v = m_grid[i];
			if (v.m_type == 'C') {
				v.m_type = 'E';
				v.m_nGroup = -1;
			}
		}

		CRange r(m_grid, fe.gv, threshold);
		while (r.next()) { 
			VOXEL& vox = r.vox();
			float dist = (r.v-fe.gv).len();
			if (dist <= threshold) {
				if (vox.m_type == 'F') { 
					minCnt++;
				} else {
					vox.m_type = 'C';
					vox.m_nGroup = -1;
				}
			}
		}
		float total = (threshold*2+1)*(threshold*2+1)*(threshold*2+1);
		int maxTotal = maxCnt(threshold);
		total = maxTotal;
		int Xsize = group(m_grid, seed, 'C');

		if (fp)
			fprintf(fp, "%d\t%d\t%d\t%d\t%f\t%d\n", threshold, threshold*2+1, minCnt, maxTotal, 1.0-minCnt/total, Xsize);

		printf("%d\t%d\t%d\t%d\t%f\t%d\n", threshold, threshold*2+1, minCnt, maxTotal, 1.0-minCnt/total, Xsize);
		return (minCnt >= m_fillCnt);
	}

	void rebuildWorld() {
		//if (m_polygons.size() == 0) return;
		m_worldBox.freeAll();
		m_worldBox.m_minSize = g_minOctreeBoxSize;
		m_worldBox.BuildWorld(m_atoms);
	}

	inline void step() { return; }
	void resize(int width, int height) { m_width = width; m_height = height; }

	inline void calculate_matrix() {
#ifdef _USE_OPENGL_
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(m_f, 1.0f, m_zNear, m_zFar);
		gluLookAt(m_vx, m_vy, m_vz, m_dx, m_dy, m_dz, m_ux, m_uy, m_uz);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(m_ox, m_oy, m_oz);
		glScalef(m_scale, m_scale, m_scale);
		glTranslatef(center.X, center.Y, center.Z);
		glRotatef(m_rx, 1.0f, 0.0f, 0.0f);
		glRotatef(m_ry, 0.0f, 1.0f, 0.0f);
		glRotatef(m_rz, 0.0f, 0.0f, 1.0f);
		glTranslatef(-center.X, -center.Y, -center.Z);
#endif
	}

	inline void drawWorld() {
#ifdef _USE_OPENGL_
		glPushMatrix();
		glDisable(GL_LIGHTING);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(m_f, 1.0f, m_zNear, m_zFar);
		gluLookAt(m_vx, m_vy, m_vz, m_dx, m_dy, m_dz, m_ux, m_uy, m_uz);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		m_worldBox.drawBox(CVector(m_x0, m_y0, m_z0), CVector(m_x1, m_y1, m_z1));
		glPopMatrix();
#endif
	}

	inline void draw() {
#ifdef _USE_OPENGL_
		int i;
		if( time(NULL) - m_lastTime >= 1 )	{
			m_lastTime = time(NULL);
			m_nFPS = m_nFrames;
			m_nFrames = 0;

			printf("FPS: %-8d #Polygons: %-8d\r", m_nFPS, g_polygonNum);
		}	m_nFrames++;

		if (m_atoms.size()==0) return;

		glColor3f(1.0f, 1.0f, 1.0f);
		if (g_viewPDB) {
			int num_polygons = (int)m_atoms.size();
			glBegin (GL_POINTS);
			for (i = 1; i < num_polygons; i++) {
				ATOM& atom = m_atoms[i];
				glVertex3f(atom.v.x, atom.v.y, atom.v.z);
			}

			ATOM& atom = m_atoms[m_seedNum];
			glColor3f(1.0f, 0, 0);
			glVertex3f(atom.v.x, atom.v.y, atom.v.z);
			glColor3f(1.0f, 1.0f, 1.0f);
			glEnd ();
		}
		
		glColor3f(1.0f, 1.0f, 1.0f);
		if (1) {
			glBegin (GL_POINTS);
			int w,d, h;
			float fw, fd, fh, unit=g_unit;
			CVector uCenter(-center.X, -center.Y, -center.Z);
			float half = m_grid.width/2;
			half*=g_unit;
			uCenter = CVector(half,half,half)+uCenter;
			i=0;
			for (fd=0,d=0; d<m_grid.depth; d++, fd+=unit) {
				for (fw=0,w=0; w<m_grid.width; w++, fw+=unit) {
					for (fh=0,h=0; h<m_grid.height; h++, fh+=unit) {
						VOXEL& vox = m_grid[i++];
						char type = vox.m_type;
						if (g_viewGrid) {
							if (type == 'F') {
								CVector v(fh, fw, fd);
								v=v-uCenter;
								glVertex3f(v.X, v.Y, v.Z);
							}
						}
						if (type == 'C') {
							CVector v(fh, fw, fd);
							v=v-uCenter;
							if (vox.m_nGroup == -1) {
								glColor3f(1.0f, 0, 0);
							} else {
								glColor3f(1.0f, 1.0f, 1.0f);
							}
							glVertex3f(v.X, v.Y, v.Z);
						}
					}
				}
			}
			glEnd ();
		}
#endif
	}

	int	m_width, m_height;
	int cmdIndex; char cmdBuffer[80], para[80];
#ifdef _USE_OPENGL_
	GLUI_String filename;
#endif
	float m_x0, m_y0, m_z0, m_x1, m_y1, m_z1, m_scale, m_ox, m_oy, m_oz, m_rx, m_ry, m_rz, m_f, m_vx, m_vy, m_vz, m_dx, m_dy, m_dz, m_ux, m_uy, m_uz, m_zNear, m_zFar;
	float m_lx, m_ly, m_lz;
	CVector minV, maxV, center;
	int windowID;
	bool m_bPaused;
	float m_ambientR, m_ambientG, m_ambientB;
	float m_diffuseR, m_diffuseG, m_diffuseB;
	char buffer[200]; 
	valarray<ATOM> m_atoms;
	CGrid m_grid;
	int	m_seedNum, m_fillCnt;

	CWorldCube m_worldBox;
	CVector m_minCorner, m_maxCorner, m_objSize;
	int		m_nFPS, m_nFrames;
	time_t	m_lastTime;
};

#endif

