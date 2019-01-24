/*--------------------------------------------------
 * main function to test reconstruction of zernike
 * discriptors. 
 * -------------------------------------------------*/

#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <GL/glut.h>

using namespace std;

#include "ZernikeDescriptor.h"
#include "Box.h"


ZernikeDescriptor<double, float>::ComplexT3D recon_result;

//#define _USE_OPENGL_
#ifdef _USE_OPENGL_

void initRendering(void)
{
	//glEnable(GL_DEPTH_TEST);

	glPointSize(5);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST );
}


float r = 0;

void display(void)
{
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0, 0, -0.2);
   	glRotatef(r, 1.0f, 2.0f, 3.0f); r+=0.1f;
	
	// clear the rendering window
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	
	
	double threshold = 0.0; 
	
	double xD = recon_result.size();
	double yD = recon_result[0].size();
	double zD = recon_result[0][0].size();

	glBegin(GL_POINTS);
	for(int z = 0; z < zD ; z++ )
	{
		for(int y = 0; y < yD; y++ )
		{
			for(int x = 0; x < xD; x++)
			{
				if( recon_result[x][y][z].real() != threshold ) 
				{
					if(recon_result[x][y][z].real() > 0)
							// set drawing color to white
						glColor3f( 1.0, 1.0, 1.0 );
					else 
						glColor3f(1.0, 0, 0);

					glVertex3f(x/(float)xD, y/(float)yD, z/(float)zD);
					
				}
			}
		}
	}
	glEnd();

	glutSwapBuffers();
	glutPostRedisplay();
}

#endif		// _USE_OPENGL_

int print3D()
{

	cout << "in print3P();\n\n";
	
	int dim = 10;
	int center = 5;
	int vol = 20;
	int min = 1;
	int max = 10;
	

	// test on sphere 
	double result[dim*dim*dim];
	
	cout << "here\n";
	for (int a = 0; a < dim; a++)
	{
		for(int b = 0; b < dim; b++)
		{	
			for(int c = 0; c < dim; c++)
			{
				if( ((a - center)*(a - center) + (b - center)*(b - center)
						+(c-center)*(c-center)) < vol)
				{
					result[(c*dim + b)*dim + a] = 1000.0;
				}else {
					result[(c*dim + b)*dim + a] = 0.0;
				}
			}
		}
	}
	
	ZernikeDescriptor<double, float> zd(result, dim, max);
	
	ZernikeMoments<float,  double> zm;
	
	//initialize the globle variable for results	
	recon_result.resize(dim);
	for(int i = 0; i < dim; i++){
		recon_result[i].resize(dim);
		for(int j = 0; j < dim; j++){
			recon_result[i][j].resize(dim);
		}
	}

	cout << "\nsize of .size() :" << recon_result.size();
	cout << "\nsize of [0].size() :" << recon_result[1].size();
	cout << "\nsize of [][].size() :" << recon_result[1][1].size() << endl;
	cout << "\nrecon function\n";
	
	zd.Reconstruct("PrintGrid.txt", recon_result, min, 100);
		
	cout << endl;
	
	// find difference by number of voxel points.
	
	zd.SaveInvariants("data.inv");

	return (dim);
}


int print3D(char* filename, int order)
{
	cout << "in print3D(" << filename << " , " << order << ")" << endl;
	
	int dim = 100;
	
	//ofstream fout("reconstruct.txt");

	// read file
	
	cout << "read file\n";
	ZernikeDescriptor<double, char> zd(filename, order);
	
	
	//initialize the globle variable for results	
	cout <<"initialize the  recon_result\n";
	recon_result.resize(dim);
	for(int i = 0; i < dim; i++){
		recon_result[i].resize(dim);
		for(int j = 0; j < dim; j++){
			recon_result[i][j].resize(dim);
		}
	}

	cout << "\nsize of .size() :" << recon_result.size();
	cout << "\nsize of [0].size() :" << recon_result[1].size();
	cout << "\nsize of [][].size() :" << recon_result[1][1].size() << endl;
	cout << "\nrecon function\n";

	char buf[100];
	char* ch_ptr = strrchr(filename, '.');
	char *reconfile;
	char *invfile;
	if(ch_ptr)
	{
		strncpy (buf, filename, ch_ptr - filename);
		buf[ch_ptr - filename] = '\0';
	}
	reconfile = buf;
	invfile = buf;
	
	strcat(reconfile, ".rc");
	
	cout << "Saving reconstructed grid file : " << reconfile << endl;
	
	zd.Reconstruct(reconfile , recon_result, 1, order, 1 ,order );
	

	cout << endl;
	
	// find difference by number of voxel points.
	strcat(invfile, ".inv");
	zd.SaveInvariants(invfile);
	

	return (dim);
}


int main(int argc, char * argv[])
{

	cout << "in main\n";
	
	if(argc == 3){
		cout << "three arguments \n";
		print3D(argv[1], atoi( argv[2] ));
	}else if(argc == 1) {
		cout << "one argument\n";
		print3D();
	}else {
		cout << "number of arguments wrong!!!:\n";
		cout << "./recon or ./recon *.inv #order\n";
	}
	

#ifdef _USE_OPENGL_
	
	int windowsize = 512;
	
	// initialize GLUT
	glutInit( &argc, argv );

	// set RGBA mode with double and depth buffers
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
	
	// set window size 
	glutInitWindowSize(windowsize, windowsize);

	//
	//glClearColor(0,0,0,1);

	glutCreateWindow(argv[0]);
	
	// 640x480, 16bit pixel depth, 60Hz refresh rate
	//glutGameModeString("640x480:16@60");

	// start fullscreen game mode
	//glutEnterGameMode();


	// Initialize OpenGL as we like it
	initRendering();
	
	// setup window callbacks
	glutDisplayFunc( display );
//	glutKeyboardFunc( keyboard );
//	glutIdleFunc( idle );


	
	//enter main loop
	glutMainLoop();

#endif // _USE_OPENGL_

	return 1;

}
