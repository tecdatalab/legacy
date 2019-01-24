// Misc. functions that have nothing to do with
// drawing a particular type of structure and can be
// used for any OpenGL application

#ifndef _OPENGL_UTILS_H_
#define _OPENGL_UTILS_H_

#include "opengl_utils.h"

#include <cstdarg>
#include <cstdio>

#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

// Naive function that returns a 3D pointer to a RGB color from a list
// defined in the .cc file. By no means, thread-safe, it's just a convenience
// when several elements need to be drawn and they just need to have
// different colors. When the list finishes it wraps around.
GLfloat* get_a_color();

void opengl_printv(va_list args, const char* format);
void opengl_print(const char* format, ...);
void opengl_printAt(int x, int y, const char* format, ...);

#endif
