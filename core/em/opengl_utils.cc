// Misc. functions that have nothing to do with
// drawing a particular type of structure and can be
// used for any OpenGL application

#include "opengl_utils.h"

// Naive way of tracking the next color to be returned by get_a_color.
// Nothing fancy is needed since it's called by a single threaded drawing
// window.
size_t current_color_index = 0;

#define TOTAL_COLORS 20

// List of colors used by get_a_color
// From http://en.wikipedia.org/wiki/List_of_colors
GLfloat color_list[TOTAL_COLORS][3] =
    {{0.4, 1, 0}, // Bright green
     {0.16, 0.32, 0.75}, // Cerulean blue
     {0.55, 0, 0.55}, // Dark magenta
     {0.89, 0.35, 0.13}, // Flame
     {0.89, 1, 0}, // Lemon lime
     {0.69, 0.77, 0.87}, // Light steel blue
     {0, 0.2, 0.13}, // Dark green
     {0.8, 0.47, 0.13}, // Ochre
     {0.5, 0, 0.5}, // Patriarch
     {0.58, 0.77, 0.45}, // Pistachio
     {0.89, 0.04, 0.36}, // Raspberry
     {0.98, 0.85, 0.37}, // Royal yellow
     {0.06, 0.32, 0.73}, // Saphire
     {0.75, 0.75, 0.75}, // Silver
     {0, 0.5, 0.5}, // Teal
     {1, 0.14, 0}, // Scarlet
     {1, 0.8, 0}, // USC Gold
     {0, 1, 1}, // Spanish sky blue
     {0.96, 0.87, 0.7}, // Wheat
     {0.91, 1, 1} // Bubble gum
    };

GLfloat* get_a_color()
{
  if(current_color_index >= TOTAL_COLORS) {
    current_color_index = 0;
  }
  return color_list[current_color_index++];
}

void opengl_printv(va_list args, const char* format)
{
	// We don't expect to print more than 8K characters
  	char buf[8192];
  	char* ch=buf;
  	vsnprintf(buf,8192,format,args);
  	while (*ch)
    	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*ch++);
}

void opengl_print(const char* format, ...)
{
  va_list args;
  va_start(args,format);
  opengl_printv(args,format);
  va_end(args);
}

void opengl_printAt(int x, int y, const char* format, ...)
{
  va_list args;
  glWindowPos2i(x,y);
  va_start(args,format);
  opengl_printv(args,format);
  va_end(args);
}

