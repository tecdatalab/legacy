#ifndef GRIDER_H
#define GRIDER_H

#pragma once
#include <stdio.h>
#include <list>
#include <vector>
#include <valarray>
#include "time.h"
#include "Grep.h"
#include <math.h>

using namespace std;

#define CSTR(str) ((char*)str.c_str())
#define VERIFY(cond, msg) { if (!(cond)) { printf("VERIFY %s FAILED @%d in %s", msg, __LINE__, __FILE__); scanf(" "); exit(-1); } }

class Vecf
{
public:
	Vecf(void) {x=y=z=0;};
	Vecf(float _x, float _y, float _z) { x=_x; y=_y; z=_z; };
	~Vecf(void) {};
	float x, y, z;
	Vecf operator +(Vecf& v) { return Vecf(v.x+x, v.y+y, v.z+z); }
	Vecf operator -(Vecf& v) { return Vecf(x-v.x, y-v.y, z-v.z);	}
	Vecf operator /(float f) { return Vecf(x/f, y/f, z/f);	}
	bool operator ==(Vecf& v) { return (v.x==x && v.y==y && v.z==z); }
	float len() { return sqrt(x*x+y*y+z*z); }
	float span() { return sqrt(x*x+y*y+z*z); }
};

struct _ATOM {
  char Atom[6]; // 1 - 6 Record name "ATOM "
  char serial[5]; // 7 - 11 Integer serial Atom serial number.
  char dummy12;
  char name[4]; // 13 - 16 Atom name Atom name.
  char altLoc; // 17   character altLoc Alternate location indicator.
  char resName[3]; // 18 - 20 Residue name resName Residue name.
  char dummy21;
  char chainID; // 22   character chainID Chain identifier.
  char resSeq[4]; // 23 - 26 Integer resSeq Residue sequence number.
  char iCode; // 27 A  char iCode Code for insertion of residues.
  char dummy28[3];
  char x[8]; // 31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms.
  char y[8]; // 39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms.
  char z[8]; // 47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms.
  char occupancy[6]; // 55 - 60 Real(6.2) occupancy Occupancy.
  char tempFactor[6]; // 61 - 66 Real(6.2) tempFactor Temperature factor.
  char dummy67[6];
  char segID[4]; // 73 - 76 LString(4) segID Segment identifier, left-justified.
  char element[2]; // 77 - 78 LString(2) element Element symbol, right-justified.
  char charge[2]; // 79 - 80 LString(2)   charge   charge on the atom.
};

struct ATOM {
  int	serial;
  char	name[5];
  char	resName[4];
  char	chainID;
  int	resSeq;
  char	iCode;
  Vecf	v;
  float radius;
  bool	bHetAtom;

  bool	bSurface;
  Vecf	gv;	//grid vector
  int	lv;	//lowest visualbility
  char	group;
  char 	binded;
};

struct VOXEL {
	int		m_atomserial, m_nGroup, m_nCubeSize;
	char	m_type;
};

class CGrid {
public:
	CGrid() { height=width=depth=0; m_size=0; m_grid=NULL; }
	~CGrid() { if (m_grid) free(m_grid); };
	void resize(int w, int h, int d) { 
		VERIFY(h==w, "h==w"); VERIFY(h==d,"h==d"); 
		m_size=w*h*d; m_grid = (VOXEL*)malloc(m_size*sizeof(VOXEL)); 
		if (!m_grid) { printf("OUTOFMEM:%d=>%d\n", w, m_size*sizeof(VOXEL)); exit(-1); }
		int i;
		for (i=0; i<m_size; i++) {
			VOXEL*p=&m_grid[i];
			p->m_atomserial=-1;
			p->m_type = 'E';
			p->m_nGroup=-1;
			p->m_nCubeSize=0;
		}
		width=w; height=h; depth=d; 
	};
	VOXEL& vox(Vecf v) {
		VERIFY(isValid((int)v.x,(int)v.y,(int)v.z), "voxV");
		return *(m_grid+((int)v.x+(int)v.y*width+(int)v.z*height*width));
	}
	VOXEL& vox(int x, int y, int z) {
		VERIFY(isValid(x,y,z), "vox");
		return *(m_grid+(x+y*width+z*height*width));
	}
	void trim(Vecf& v) {
		while ((int)v.x<0) v.x++; while ((int)v.y<0) v.y++; while ((int)v.z<0) v.z++;
		while ((int)v.x>=width) v.x--; while ((int)v.y>=height) v.y--; while ((int)v.z>=depth) v.z--;
		VERIFY(isValid((int)v.x,(int)v.y,(int)v.z), "trim");
	}
	inline bool isValid(int x, int y, int z) { return ((x>=0) && (y>=0) && (z>=0) && (x<width) && (y<height) && (z<depth)); }
	inline bool isValid(Vecf v) { return isValid((int)v.x,(int)v.y,(int)v.z); }

	inline VOXEL& operator[] (const int index) { return *(m_grid+index); }
	int height, width, depth;
private:
	VOXEL* m_grid;
	int m_size;
};

class CRange {
public:
	CRange(CGrid& grid) : m_grid(grid) { 
		Vecf start = Vecf(0, 0, 0);
		Vecf end = Vecf((float)m_grid.width-1, (float)m_grid.height-1, (float)m_grid.depth-1);
		m_center = Vecf(-999,-999,-999); init(start, end);
	}
	CRange(CGrid& grid, Vecf center, float radius): m_grid(grid) { 
		VERIFY(radius>0, "radius>0");
		VERIFY(m_grid.isValid(center), "center");
		Vecf start = Vecf(center.x-radius,center.y-radius,center.z-radius); grid.trim(start);
		Vecf end = Vecf(center.x+radius,center.y+radius,center.z+radius); grid.trim(end);
		m_center = center;	init(start, end);
	}
	bool next() { 
		if (!m_hasNext) 
			return false; 
		if (!m_bFirst) 
			++(*this); 
		m_bFirst=false; 
		VERIFY(!(v==m_center), "next");
		return m_hasNext; 
	};
	void operator++() {
		do {
			do { //VERIFY(!m_hasNext, "m_hasNext");
				if (++v.z > m_end.z) v.z = m_start.z; else break;
				if (++v.y > m_end.y) v.y = m_start.y; else break;
				if (++v.x > m_end.x) { m_hasNext = false; break; }
			} while (!m_grid.isValid(v));
		} while (v==m_center);
	};
	void rewind() {	init(m_start, m_end); }

	VOXEL& vox() { return m_grid.vox(v); }
	Vecf v;
private:
	void init(Vecf start, Vecf end) {
		m_bFirst = true;
		m_hasNext = true;	
		m_start = start;
		m_end = end;
		v = m_start;
		while ((!m_grid.isValid(v) && m_hasNext) || (v==m_center) )
			++(*this);
	}
	CGrid& m_grid;
	Vecf m_start, m_end, m_center;
	bool m_bFirst, m_hasNext;
	VOXEL dummy;
};

class CReader
{
protected:
	valarray<string> m_lines;
public:
	CReader(void) {};
	~CReader(void) {};
	int read(char* filename) {
		char buf[100]; size_t line;
		list<string> lst;
		//printf("Reading %s ", filename);
		FILE* fp = fopen(filename, "r"); if (!fp) return 0;
		while (fgets(buf, 100, fp) != NULL) {
			lst.push_back(buf);
		}
		fclose(fp);
		m_lines.resize(lst.size());
		list<string>::iterator iter = lst.begin();
		for (line=0; line<lst.size(); line++, iter++) {
			m_lines[line] = (*iter).c_str();
		};
		return 1;
	};
	int size() { return (int)m_lines.size(); }
	const char* operator[](int index) {
		return m_lines[index].c_str();
	}
	int exist(char* filename) {
		FILE* fp = fopen(filename, "r");
		if (!fp) return 0;
		fclose(fp); return 1;
	}
};

int group(CGrid& grid, Vecf seed, char type) {
	int count = 0;
	list<Vecf> neighbours; 
	VOXEL& vox=grid.vox(seed);
	if (vox.m_nGroup == -1) {
		vox.m_nGroup = 1;
		neighbours.push_front(seed);
		while(1) {
			if ( neighbours.size() == 0 ) break;
			Vecf curVec = neighbours.front(); neighbours.pop_front();
			CRange r(grid, curVec, 2);
			while (r.next()) {
				VOXEL& vox2 = r.vox();
				if ((vox2.m_type == type) && vox2.m_nGroup == -1) {
					count++;
					vox2.m_nGroup = vox.m_nGroup;
					neighbours.push_front(r.v);
				}
			}
		}
	}

	return count;
}

#endif
