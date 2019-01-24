#include "pdb_painter.h"

pdb_painter:: pdb_painter(pdb* in_protein) : protein(in_protein) {
}

void pdb_painter::draw_calpha()
{
	glPushMatrix();
  glLineWidth(5);
  glBegin(GL_LINE_STRIP);

  for (size_t atom_index = 0; atom_index < protein->atoms.size();
       atom_index++) {
    atom current_atom = (protein->atoms)[atom_index];
    if (!current_atom.atype.compare("CA")) {
      // Assign any normal just to have some shadow. We don't really worry
      // too much about the direction at this point.
      glNormal3f(1.0, 0.0, 0.0);
      glVertex3f(current_atom.axyz[0], current_atom.axyz[1],
                 current_atom.axyz[2]);
    } 
  }
  glEnd();
	glPopMatrix();
}
