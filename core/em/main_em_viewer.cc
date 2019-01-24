#include "correlation_results.h"
#include "density_lattice.h"
#include "main_em_correlation_options.h"
#include "marching_cubes_surface.h"
#include "mrc.h"
#include "opengl_utils.h"

#include "../proteinstructure/pdb.h"
#include "../proteinstructure/pdb_painter.h"
#include "../proteinstructure/rmsd.h"
#include "../proteinstructure/transformable_pdb.h"

/*
*/
#include <cmath>
#include <iostream>
#include <vector>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>		   // Open Graphics Library (OpenGL) header
#include <GL/glut.h>	   // The GL Utility Toolkit (GLUT) Header
#endif

#define KEY_ESCAPE 27

typedef struct {
    int width;
	int height;
	char* title;

	float field_of_view_angle;
	float z_near;
	float z_far;
    // Look at specification
    GLdouble eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz;
    // Eye direction
} glutWindow;
glutWindow win;

correlation_results results;
mrc_reader* ref_reader;
mrc_reader* fitted_reader;
density_lattice* ref_lattice;
density_lattice* fitted_lattice;
marching_cubes_surface* ref_surface;
marching_cubes_surface* fitted_surface;
double isovalue = 0;
size_t results_index = 0;
double current_rmsd = 0.0f;
transformable_pdb ref_pdb;
transformable_pdb fitted_pdb;
pdb_painter* ref_pdb_painter;

density_vector min_bound, max_bound;
density_vector fitted_min_bound, fitted_max_bound;
//unsigned long max_overlap;

int azimuth = 90;
int elevation = 90;

float aspect = 1;
float dim = 200;

// Rotation with the arrows, in degrees
#define ROTATION 5
#define degree_cos(th) cos(M_PI/180*(th))
#define degree_sin(th) sin(M_PI/180*(th))

#define AXIS_LINE_LENGTH 25
bool enable_axes = 1;

// This flag determines if separate surfaces are displayed for each
// prediction (true) or if a point at the center of each loaded prediction
// is shown
bool show_separate_surfaces = 1;

bool show_reference_map = 1;

bool show_unmodified_pdbs = 1;

void project() 
{
	glMatrixMode(GL_PROJECTION);
  	glLoadIdentity();

    gluPerspective(win.field_of_view_angle,aspect,dim/4,4*dim);

 	glMatrixMode(GL_MODELVIEW);
  	glLoadIdentity();
}
/*
 *  drawLight 
 *  ------
 *  Draws the light
 */
#define DEF_DISTANCE 50
#define DEF_AMBIENT 45
#define DEF_DIFFUSE 100
#define DEF_SPECULAR 10
#define DEF_L_Y 0
#define DEF_L_PH 0

int light_distance=DEF_DISTANCE;
int ambient=DEF_AMBIENT;
int diffuse=DEF_DIFFUSE;
int specular=DEF_SPECULAR;
float lightY=DEF_L_Y;
int lightPh=DEF_L_PH;

void drawLight()
{
	/*  Translate intensity to color vectors */
    GLfloat Ambient[]   = {0.01*ambient ,0.01*ambient ,0.01*ambient ,1.0};
    GLfloat Diffuse[]   = {0.01*diffuse ,0.01*diffuse ,0.01*diffuse ,1.0};
    GLfloat Specular[]  = {0.01*specular,0.01*specular,0.01*specular,1.0};
    GLfloat Position[]  = {light_distance*degree_sin(lightPh),lightY,
							light_distance*degree_cos(lightPh),1.0};

    /*  Set ambient, diffuse, specular components and position of light 0 */
    /*
      Light uses the Phong model
      Once light is enabled, colors assigned by glColor* isn't used
      Ambient is light that's been scattered by environment that its direction is impossible to determine
      Diffuse is is light that comes from one direction, so it's brighter if it comes squarely on surface rather than glances off
      Specular is light that comes from a particular direction and bounces off in preferred direction
      Position is the position of our light. In this case it is the same as the sphere.
     */
    glLightfv(GL_LIGHT0,GL_AMBIENT, Ambient);
    glLightfv(GL_LIGHT0,GL_DIFFUSE, Diffuse);
    glLightfv(GL_LIGHT0,GL_POSITION,Position);

    glLightfv(GL_LIGHT0,GL_SPECULAR,Specular);
    /*  Use glColorMaterial when you need to change a single material parameter for most vertices
	in your scene */
    /*  glColorMaterial sets ambient and diffuse color materials */
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    /*  Now glColor* changes ambient and diffuse reflection */
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING); 
    glEnable(GL_LIGHT0);
}

void drawAxes()
{
	if (enable_axes) {
      glLineWidth(2.5);
	    glDisable(GL_LIGHTING);
    	glColor3f(0.0f, 0.0f , 0.0f);
	    glBegin(GL_LINES);
	    glVertex3d(0,0,0);
	    glVertex3d(AXIS_LINE_LENGTH,0,0);
	    glVertex3d(0,0,0);
	    glVertex3d(0,AXIS_LINE_LENGTH,0);
	    glVertex3d(0,0,0);
	    glVertex3d(0,0,AXIS_LINE_LENGTH);
	    glEnd();
		/*  Label axes */
    	glRasterPos3d(AXIS_LINE_LENGTH,0,0);
	    opengl_print("X");
    	glRasterPos3d(0,AXIS_LINE_LENGTH,0);
	    opengl_print("Y");
    	glRasterPos3d(0,0,AXIS_LINE_LENGTH);
	    opengl_print("Z");
    	glEnable(GL_LIGHTING);
	}
}

// For results[index] it draws either a full surface (if show_surface=true)
// or a point centered
void draw_fitted_structure(size_t index, bool show_surface)
{
	glPushMatrix();

	glTranslatef((results[index].x + ref_lattice->origin_x())
                   * ref_lattice->voxel_length_x(),
					     (results[index].y + ref_lattice->origin_y())
                   * ref_lattice->voxel_length_y(),
					     (results[index].z + ref_lattice->origin_z())
                   * ref_lattice->voxel_length_z());
	glRotatef(results[index].alpha_z, 0, 0, 1);
	glRotatef(results[index].beta_y, 0, 1, 0);
	glRotatef(results[index].gamma_x, 1, 0, 0);

	if(show_surface) {
		glColor4f(1,1,0,0.5);
		fitted_surface->draw(isovalue);
	  glPopMatrix();
    // Draw the fitted PDB by manually applying the transformation.
    // This is analogous to using the OpenGL transformation but, for
    // testing purposes, it's better to use the manually transformed version.
    glColor4f(0.16, 0.32, 0.75, 0);
    transformable_pdb pdb_copy = fitted_pdb;

    pdb_copy.translate(
      -fitted_reader->get_x_origin() * fitted_lattice->voxel_length_x(),
      -fitted_reader->get_y_origin() * fitted_lattice->voxel_length_y(),
      -fitted_reader->get_z_origin() * fitted_lattice->voxel_length_z());

    point_transformation transformation(results[index].alpha_z,
        results[index].beta_y, results[index].gamma_x,
        (results[index].x + ref_lattice->origin_x())
            * ref_lattice->voxel_length_x(),
        (results[index].y + ref_lattice->origin_y())
            * ref_lattice->voxel_length_y(),
        (ref_lattice->origin_z() + results[index].z)
            * ref_lattice->voxel_length_z());
    pdb_copy.apply_point_transformation(&transformation);

    pdb_painter fitted_pdb_painter(&pdb_copy);
    fitted_pdb_painter.draw_calpha();

    // Update the RMSD between the fitted PDB and the reference
    transformable_pdb fitted_matches;
    transformable_pdb ref_matches;

    pdb_copy.get_matching_calphas(&ref_pdb, true, &(fitted_matches.atoms),
        &(ref_matches.atoms));
/*
    fitted_matches.scale(fitted_lattice->voxel_length_x(),
                         fitted_lattice->voxel_length_y(),
                         fitted_lattice->voxel_length_z());

    ref_matches.scale(ref_lattice->voxel_length_x(),
                      ref_lattice->voxel_length_y(),
                      ref_lattice->voxel_length_z());
  */  
    current_rmsd = calculate_allatom_rmsd(fitted_matches.atoms,
                                          ref_matches.atoms);
    // Test if the min and max bounds after transformation are correct.
    // Just draw a diagonal across the structure.
    density_vector transformed_boundary_min, transformed_boundary_max;
    point_transformation t1(0, 0, 0, 
      -fitted_reader->get_x_origin() * fitted_lattice->voxel_length_x(),
      -fitted_reader->get_y_origin() * fitted_lattice->voxel_length_y(),
      -fitted_reader->get_z_origin() * fitted_lattice->voxel_length_z());
    point_transformation_sequence t;
    t.add(t1);
    t.add(transformation);

    fitted_lattice->get_boundaries_after_transformation(t, fitted_min_bound,
        fitted_max_bound, &transformed_boundary_min, &transformed_boundary_max);

    unsigned long vx_min, vy_min, vz_min, vx_max, vy_max, vz_max;
    ref_lattice->real_coordinate_to_voxel(transformed_boundary_min.x,
        transformed_boundary_min.y, transformed_boundary_min.z, &vx_min,
        &vy_min, &vz_min);
    ref_lattice->real_coordinate_to_voxel(transformed_boundary_max.x,
        transformed_boundary_max.y, transformed_boundary_max.z, &vx_max,
        &vy_max, &vz_max);

cout << "Voxels: " << vx_min << " " << vy_min << " " << vz_min
     << " - " << vx_max << " " << vy_max << " " << vz_max << "\n";
double test_x, test_y, test_z;
ref_lattice->voxel_to_real_coordinate(vx_min, vy_min, vz_min, &test_x, &test_y, &test_z);
ref_lattice->voxel_to_real_coordinate(vx_max, vy_max, vz_max, &test_x, &test_y, &test_z);

    glLineWidth(5);
    glColor3f(0.0f, 0.0f , 0.0f);
	  glBegin(GL_LINES);
	  glVertex3f(transformed_boundary_min.x, transformed_boundary_min.y,
               transformed_boundary_min.z);
	  glVertex3f(transformed_boundary_min.x, transformed_boundary_min.y,
               transformed_boundary_max.z);
	  glVertex3f(transformed_boundary_max.x, transformed_boundary_min.y,
               transformed_boundary_max.z);
	  glVertex3f(transformed_boundary_max.x, transformed_boundary_min.y,
               transformed_boundary_min.z);
	  glVertex3f(transformed_boundary_min.x, transformed_boundary_min.y,
               transformed_boundary_min.z);
	  glEnd();
	  glBegin(GL_LINES);
	  glVertex3f(transformed_boundary_min.x, transformed_boundary_max.y,
               transformed_boundary_min.z);
	  glVertex3f(transformed_boundary_min.x, transformed_boundary_max.y,
               transformed_boundary_max.z);
	  glVertex3f(transformed_boundary_max.x, transformed_boundary_max.y,
               transformed_boundary_max.z);
	  glVertex3f(transformed_boundary_max.x, transformed_boundary_max.y,
               transformed_boundary_min.z);
	  glVertex3f(transformed_boundary_min.x, transformed_boundary_max.y,
               transformed_boundary_min.z);
	  glEnd();

    double test_overlap;
    double test_correlation = ref_lattice->calculate_correlation(
        *fitted_lattice, t, transformed_boundary_min, transformed_boundary_max,
        isovalue, &test_overlap);

    cout << "Correlation " << test_correlation
         << " Overlap " << test_overlap << endl;
    
	} else {
		glColor4f(0,1,1,0.5);
		glPointSize(8.0);
    glTranslatef(
        -fitted_reader->get_x_origin() * fitted_lattice->voxel_length_x(),
        -fitted_reader->get_y_origin() * fitted_lattice->voxel_length_y(),
        -fitted_reader->get_z_origin() * fitted_lattice->voxel_length_z());
    glBegin(GL_POINTS);
		glVertex3f(
			((fitted_max_bound.x - fitted_min_bound.x) / 2.0f) + fitted_min_bound.x,
			((fitted_max_bound.y - fitted_min_bound.y) / 2.0f) + fitted_min_bound.y,
			((fitted_max_bound.z - fitted_min_bound.z) / 2.0f) + fitted_min_bound.z);
		glEnd();
	  glPopMatrix();
	}
}

void display() 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		     // Clear Screen and Depth Buffer
	/*  Enable Z-buffering in OpenGL */
  glEnable(GL_DEPTH_TEST);
//  glEnable (GL_BLEND);
//  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLoadIdentity();

	double eyex = -2*dim*degree_sin(azimuth)*degree_cos(elevation);
    double eyey = 2*dim*degree_sin(elevation);
    double eyez = 2*dim*degree_cos(azimuth)*degree_cos(elevation);
    /* camera/eye position, aim of camera lens, up-vector */
    gluLookAt(eyex,eyey,eyez , 0,0,0 , 0,degree_cos(elevation),0);
	
	drawAxes();
	drawLight();
	
	glPushMatrix();										  // Push the current matrix stack
/*
	glTranslatef(
		((int)(max_bound.x - min_bound.x)) / -2 - ((int)min_bound.x),
		((int)(max_bound.y - min_bound.y)) / -2 - ((int)min_bound.y),
		((int)(max_bound.z - min_bound.z)) / -2 - ((int)min_bound.z));
*/
	glColor4f(1,0,0,0.5);
    //drawDensityMap(ref_lattice);
  if(show_reference_map) {
    glPushMatrix();
	  glTranslatef(ref_lattice->origin_x() * ref_lattice->voxel_length_x(),
					      ref_lattice->origin_y() * ref_lattice->voxel_length_y(),
					      ref_lattice->origin_z() * ref_lattice->voxel_length_z());
  	ref_surface->draw(isovalue);
    glPopMatrix();
    glColor4f(0.5, 0, 0.5, 0);
    ref_pdb_painter->draw_calpha();
  }

  if(show_unmodified_pdbs) {
    //draw_pdbs();
  }

	//cerr << "Fitted: " << results[0].x << " " << " " << results[0].y << " " << results[0].z << endl;
	if(show_separate_surfaces) {
		draw_fitted_structure(results_index, true);
	} else {
		for(size_t index = 0; index < results.size(); index++) {
			draw_fitted_structure(index, false);
		}
	}

	glPopMatrix();										  // Pop the current matrix stack
	glColor4f(0,0,0,0.5);
	opengl_printAt(5,10, "X: %d Y: %d Z: %d rX: %f rY: %f rZ: %f",
					(int)results[results_index].x,
					(int)results[results_index].y,
					(int)results[results_index].z,
					(float)results[results_index].gamma_x,
					(float)results[results_index].beta_y,
					(float)results[results_index].alpha_z);
	opengl_printAt(5,30, "CC: %f Overlap: %f RMSD: %f, Prediction# %d/%d",
					(float)results[results_index].correlation,
					(float)results[results_index].overlap,
          (float)current_rmsd,
					results_index + 1, results.size());

	glFlush();
	glutSwapBuffers();
}

void reshape(int width,int height)
{
  	aspect = (height>0) ? (double)width/height:1;
  	glViewport(0,0,width,height);
	project();
}

void keyboard ( unsigned char key, int mousePositionX, int mousePositionY )		
{ 
  switch ( key ) 
  {
    case KEY_ESCAPE:        
      exit(0);   
      break;
    case 'n':
      results_index++;
      if (results_index >= results.size()) {
        results_index = 0;
      }
      break;
		case 'p':
      if (results_index > 0) {
        results_index--;
      } else {
        results_index = results.size() - 1;
      }
      break;
    case 'w':
      ref_surface->set_wireframe(!ref_surface->is_wireframe_only());
      break;
    case 'm':
      show_separate_surfaces = !show_separate_surfaces;
      break;
    case 'r':
      show_reference_map = !show_reference_map;
      break;
    case 'f':
      show_unmodified_pdbs = !show_unmodified_pdbs;
    default:      
    		break;
  }
  glutPostRedisplay();
}

void specialKeyboard(int key, int x, int y)
{
    int mod = glutGetModifiers();
	switch (key)
	{
		case GLUT_KEY_RIGHT:
			azimuth += ROTATION;
			break;
		case GLUT_KEY_LEFT:
			azimuth -= ROTATION;
			break;
		case GLUT_KEY_UP:
			if (mod == GLUT_ACTIVE_SHIFT) { //zoom
				win.field_of_view_angle--;
			} else {
				elevation += ROTATION;
			}
			break;
		case GLUT_KEY_DOWN:
			if (mod == GLUT_ACTIVE_SHIFT) { // zoom
				win.field_of_view_angle++;
			} else {
				elevation -= ROTATION;
			}
			break;
	}
  	/* If the angles circle around */
	azimuth %= 360;
	elevation %= 360;

	project();
  glutPostRedisplay();
}

int main(int argc, char **argv) 
{
  if(argc < 7) {
    cerr << "Usage: " << argv[0] << " reference_mrc_file fitted_mrc_file "
	       << "transformations_file isovalue reference_pdb fitted_pdb" << endl;
    return 0;
  }
  ref_reader = new mrc_reader(argv[1]);
  fitted_reader = new mrc_reader(argv[2]);
	results.read(argv[3]);
	isovalue = atof(argv[4]);

  ref_lattice = ref_reader->get_densities();
  fitted_lattice = fitted_reader->get_densities();

	ref_surface = new marching_cubes_surface(ref_lattice, true);
	fitted_surface = new marching_cubes_surface(fitted_lattice, true);

  read_protein(argv[5], ref_pdb);
  read_protein(argv[6], fitted_pdb);

/*  ref_pdb.scale(1.0f / ref_lattice->voxel_length_x(),
                1.0f / ref_lattice->voxel_length_y(),
                1.0f / ref_lattice->voxel_length_z());*/
/*
  ref_pdb.translate(-ref_reader->get_x_origin() * ref_lattice->voxel_length_x(),
      -ref_reader->get_y_origin() * ref_lattice->voxel_length_y(),
      -ref_reader->get_z_origin() * ref_lattice->voxel_length_z());

  fitted_pdb.translate(
      -fitted_reader->get_x_origin() * fitted_lattice->voxel_length_x(),
      -fitted_reader->get_y_origin() * fitted_lattice->voxel_length_y(),
      -fitted_reader->get_z_origin() * fitted_lattice->voxel_length_z());

  fitted_pdb.scale(1.0f / fitted_lattice->voxel_length_x(),
                   1.0f / fitted_lattice->voxel_length_y(),
                   1.0f / fitted_lattice->voxel_length_z());*/
  
  ref_pdb_painter = new pdb_painter(&ref_pdb);

  // TODO: this was changed to real values, fix it
	ref_lattice->get_boundaries_for_threshold(isovalue, &min_bound, &max_bound);
	fitted_lattice->get_boundaries_for_threshold(isovalue,
										&fitted_min_bound, &fitted_max_bound);

	//max_overlap = fitted_lattice->count_voxels_over_threshold(TMP_ISOVALUE);

	// Get high correlation translations
	// TODO: This is useless now because the idea is not to compute these
	// but to load them from a file
	//results = ref_lattice->get_significant_correlations(fitted_lattice,
	//			5, 90, TMP_ISOVALUE, 0.5, 0.5);

	// set window values
	win.width = 900;
	win.height = 900;
	win.title = (char*)"Em Viewer";
	win.field_of_view_angle = 45;
	win.z_near = 1.0f;
	win.z_far = 500.0f;
  // look at parameters
  win.eyex = 0;
  win.eyey = 0;
  win.eyez = 200;
  win.centerx = 0;
  win.centery = 0;
  win.centerz = 0;
  win.upx = 0;
  win.upy = 1;
  win.upz = 0;

	// initialize and run program
	glutInit(&argc, argv);                                      // GLUT initialization
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );  // Display Mode
	glutInitWindowSize(win.width,win.height);					// set window size
	glutCreateWindow(win.title);								// create Window
  glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glutDisplayFunc(display);									// register Display Function
	glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);								// register Keyboard Handler
    glutSpecialFunc(specialKeyboard);
	glutMainLoop();												// run GLUT mainloop


  delete ref_pdb_painter;
	delete fitted_surface;
	delete ref_surface;
	delete fitted_lattice;
  delete ref_lattice;
  delete fitted_reader;
	delete ref_reader;
	return 0;
}
