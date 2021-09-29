// cd Desktop/MTECH/Graphics/Assignments/Assignment2/code/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <glm/vec3.hpp>
#include <chrono>

#include "file_utils.h"
#include "math_utils.h"
#include "POTReader.h"
#include "texture_3d.h"

using namespace std;
using namespace std::chrono;

/********************************************************************/
/*   Variables */

const char * theProgramTitle = "Assignment2-Nikita";
int theWindowWidth = 1200, theWindowHeight = 700;
int theWindowPositionX = 40, theWindowPositionY = 40;
bool isFullScreen = false;
bool isAnimating = true;
float rotation = 0.0f;
GLuint VBO, VAO, nbo, tbo, ibo;
int np,n;
GLint gWorldLocation,gllightColor,globjectColor,gllightPos,glviewPos,gfieldvmin,gfieldvmax;
GLint ga,gb,gc,gd,side,gslice;
int setside=1;
int enable_slice=2;
ScalarField *field;

/* Constants */
const int ANIMATION_DELAY = 20; /* milliseconds between rendering */
const char* pVSFileName = "shader.vs";
const char* pFSFileName = "shader.fs";

/*Texture*/
Texture *Texture_3d=NULL,*Texture_1d=NULL;
GLuint gSampler,gSampler1,gIsovalue;
float isovalue=0.5;
float *colorTexture;


/********************************************************************
  Utility functions
 */
/*********************************************************************
/* post: compute frames per second and display in window's title bar */
void computeFPS() {
	static int frameCount = 0;
	static int lastFrameTime = 0;
	static char * title = NULL;
	int currentTime;

	if (!title)
		title = (char*) malloc((strlen(theProgramTitle) + 20) * sizeof (char));
	frameCount++;
	currentTime = glutGet((GLenum) (GLUT_ELAPSED_TIME));
	if (currentTime - lastFrameTime > 1000) {
		sprintf(title, "%s [ FPS: %4.2f ]",
			theProgramTitle,
			frameCount * 1000.0 / (currentTime - lastFrameTime));
		glutSetWindowTitle(title);
		lastFrameTime = currentTime;
		frameCount = 0;
	}
}

/****************************************************************/
static float normaliseVertexAttribute(float x, float min, float max) {
	return (2*((x-min)/(max-min))-1.0);
} 

static float normaliseScalarField(float s, float smin, float smax ) {
	return (s-smin)/(smax-smin);
}

static Vector3f getColor(float alpha) {
	Vector3f red=Vector3f(1.0f,0.0f,0.0f);
	Vector3f blue=Vector3f(0.0f,0.0f,1.0f);

	return red + (blue-red)*((alpha-field->vmin)/(field->vmax-field->vmin));
}
/*****************************************************************/

void create1dTex() {
	float alpha;
	colorTexture = new float[99];
	for(int i=0;i<99;i+=3){
		alpha=((field->vmax-field->vmin)/99)*i+field->vmin;
		Vector3f c=getColor(alpha);
		colorTexture[i]=c.x;
		colorTexture[i+1]=c.y;
		colorTexture[i+2]=c.z;
	}	
	
}

static void CreateTexture() {
	create1dTex();
	Texture_3d= new Texture(GL_TEXTURE_3D);
	Texture_3d->Load3d(field);
	Texture_1d= new Texture(GL_TEXTURE_1D);
	Texture_1d->Load1d(colorTexture);
}

static void CreateVertexBuffer() {
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	const char* OffFile = "1grm.off";
	int nv,noEdges,i,j;
	char type[3];
	FILE* input;

	float xmin,ymin,zmin,xmax,ymax,zmax;

	input = fopen(OffFile, "r");
	fscanf(input, "%s", type);
	/* First line should be OFF */
	if(strcmp(type,"OFF")) {
		printf("Not a OFF file");
		exit(1);
	}
	/* Read the no. of vertices, faces and edges */
	fscanf(input, "%d", &nv);
	fscanf(input, "%d", &np);
	fscanf(input, "%d", &noEdges);

	Vector3f *vertices=(Vector3f *)malloc(sizeof(Vector3f)*nv);
	float x,y,z;
	/* Read the vertices' location*/	
	for(i = 0;i < nv;i++) {
		fscanf(input, "%f %f %f", &x,&y,&z);
		vertices[i].x = x;
		vertices[i].y = y;
		vertices[i].z = z;
		//cout << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << "\n";
		if (i==0){
			xmin=xmax=x;
			ymin=ymax=y;
			zmin=zmax=z;
		}
		else {
			xmin=min(xmin,x);
			xmax=max(xmax,x);
			ymin=min(ymin,y);
			ymax=max(ymax,y);
			zmin=min(zmin,z);
			zmax=max(zmax,z);
		}
	}

	/*Scalar Field*/
	field = loadField("1grm.pot");

	Vector3f *Texcoor=(Vector3f *)malloc(sizeof(Vector3f)*nv);
	//Texture coordinate
	for ( i=0;i<nv;i++ ) {
		Texcoor[i].x=(vertices[i].x-field->origin[0])/(field->step[0]*(field->dim[0]-1));
		Texcoor[i].y=(vertices[i].y-field->origin[1])/(field->step[1]*(field->dim[1]-1));
		Texcoor[i].z=(vertices[i].z-field->origin[2])/(field->step[2]*(field->dim[2]-1));
	}


	/*Normalize vertices*/
	/*for(i = 0;i < nv;i++) {
		vertices[i].x = normaliseVertexAttribute(vertices[i].x,xmin,xmax);
		vertices[i].y = normaliseVertexAttribute(vertices[i].y,ymin,ymax);
		vertices[i].z = normaliseVertexAttribute(vertices[i].z,zmin,zmax);
	}*/

	for(i = 0;i < nv;i++) {
		vertices[i].x = normaliseVertexAttribute(vertices[i].x,field->origin[0],field->origin[0]+(field->dim[0]-1)*field->step[0]);
		vertices[i].y = normaliseVertexAttribute(vertices[i].y,field->origin[1],field->origin[1]+(field->dim[1]-1)*field->step[1]);
		vertices[i].z = normaliseVertexAttribute(vertices[i].z,field->origin[2],field->origin[2]+(field->dim[2]-1)*field->step[2]);
	}

	std::vector<unsigned int> faces;
	unsigned int v;
	/* Read the Polygons */	
	for(i = 0;i < np;i++) {
		/* No. of sides of the polygon (Eg. 3 => a triangle) */
		fscanf(input, "%d", &n);
		//cout << n << " ";
		
		/* read the vertices that make up the polygon */
		for(j = 0;j < n;j++) {
			fscanf(input, "%d", &v);
			//cout << v << " ";
			faces.push_back(v);
		}
		//cout <<"\n";
	}

	/*Polygon indices*/
	unsigned int *indices=(unsigned int*)malloc(sizeof(unsigned int)*np*3);
	j=0;

	for( std::vector<unsigned int>::iterator it=faces.begin(); it!=faces.end(); it++ ){
		indices[j]=*it;
		j++;
	}

	/* Calculate normals*/
	Vector3f *normals=(Vector3f *)malloc(sizeof(Vector3f)*nv);

	for (i=0;i<nv;i++) {
		normals[i]=Vector3f(0.0f,0.0f,0.0f);
	}

	int limit=np*3;
	for (i=0;i<limit;i+=3) {
		int idx0=indices[i];
		int idx1=indices[i+1];
		int idx2=indices[i+2];
		Vector3f v1 = vertices[idx1];
		v1-=vertices[idx0];
		Vector3f v2 = vertices[idx2];
		v2-=vertices[idx0];
		Vector3f normal = v1.Cross(v2);
		normal.Normalize();

		normals[idx0]+=normal;
		normals[idx1]+=normal;
		normals[idx2]+=normal;
	}

	for (i=0;i<nv;i++) {
		normals[i].Normalize();
	}

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vector3f)*nv, vertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glGenBuffers(1,&nbo);
	glBindBuffer(GL_ARRAY_BUFFER, nbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vector3f)*nv, normals, GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,0);

	
	glGenBuffers(1,&tbo);
	glBindBuffer(GL_ARRAY_BUFFER,tbo);
	glBufferData(GL_ARRAY_BUFFER,sizeof(Vector3f)*nv,Texcoor,GL_STATIC_DRAW);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,0,0);

	
	glGenBuffers(1,&ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*np*sizeof(unsigned int), indices, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	CreateTexture();
}

static void AddShader(GLuint ShaderProgram, const char* pShaderText, GLenum ShaderType) {
	GLuint ShaderObj = glCreateShader(ShaderType);

	if (ShaderObj == 0) {
		fprintf(stderr, "Error creating shader type %d\n", ShaderType);
		exit(0);
	}

	const GLchar * p[1];
	p[0] = pShaderText;
	GLint Lengths[1];
	Lengths[0] = strlen(pShaderText);
	glShaderSource(ShaderObj, 1, p, Lengths);
	glCompileShader(ShaderObj);
	GLint success;
	glGetShaderiv(ShaderObj, GL_COMPILE_STATUS, &success);
	if (!success) {
		GLchar InfoLog[1024];
		glGetShaderInfoLog(ShaderObj, 1024, NULL, InfoLog);
		fprintf(stderr, "Error compiling shader type %d: '%s'\n", ShaderType, InfoLog);
		exit(1);
	}

	glAttachShader(ShaderProgram, ShaderObj);
}


static void CompileShaders() {
	GLuint ShaderProgram = glCreateProgram();

	if (ShaderProgram == 0) {
		fprintf(stderr, "Error creating shader program\n");
		exit(1);
	}

	string vs, fs;

	if (!ReadFile(pVSFileName, vs)) {
		exit(1);
	}

	if (!ReadFile(pFSFileName, fs)) {
		exit(1);
	}

	AddShader(ShaderProgram, vs.c_str(), GL_VERTEX_SHADER);
	AddShader(ShaderProgram, fs.c_str(), GL_FRAGMENT_SHADER);

	GLint Success = 0;
	GLchar ErrorLog[1024] = {0};

	glLinkProgram(ShaderProgram);
	glGetProgramiv(ShaderProgram, GL_LINK_STATUS, &Success);
	if (Success == 0) {
		glGetProgramInfoLog(ShaderProgram, sizeof (ErrorLog), NULL, ErrorLog);
		fprintf(stderr, "Error linking shader program: '%s'\n", ErrorLog);
		exit(1);
	}

/*	glValidateProgram(ShaderProgram);
	glGetProgramiv(ShaderProgram, GL_VALIDATE_STATUS, &Success);
	if (!Success) {
		glGetProgramInfoLog(ShaderProgram, sizeof (ErrorLog), NULL, ErrorLog);
		fprintf(stderr, "Invalid shader program: '%s'\n", ErrorLog);
		exit(1);
	}*/

	glUseProgram(ShaderProgram);
	gWorldLocation = glGetUniformLocation(ShaderProgram, "gWorld");
	gllightColor = glGetUniformLocation(ShaderProgram,"lightColor");
	//globjectColor = glGetUniformLocation(ShaderProgram,"objectColor");
	gllightPos = glGetUniformLocation(ShaderProgram,"lightPos");
	glviewPos = glGetUniformLocation(ShaderProgram,"viewPos");
	gSampler = glGetUniformLocation(ShaderProgram,"gampler");
	gSampler1 = glGetUniformLocation(ShaderProgram,"gSampler1");
	gIsovalue = glGetUniformLocation(ShaderProgram,"gIsovalue");
	gfieldvmin = glGetUniformLocation(ShaderProgram,"fieldvmin");
	gfieldvmax = glGetUniformLocation(ShaderProgram,"fieldvmax");
	ga = glGetUniformLocation(ShaderProgram,"ga");
	gb = glGetUniformLocation(ShaderProgram,"gb");
	gc = glGetUniformLocation(ShaderProgram,"gc");
	gd = glGetUniformLocation(ShaderProgram,"gd");
	side = glGetUniformLocation(ShaderProgram,"side");
	gslice = glGetUniformLocation(ShaderProgram,"gslice");
}

/********************************************************************
 Callback Functions
 These functions are registered with the glut window and called 
 when certain events occur.
 */

void onInit(int argc, char * argv[])
/* pre:  glut window has been initialized
   post: model has been initialized */ {
	/* by default the back ground color is black */
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	CreateVertexBuffer();
	CompileShaders();

	/* set to draw in window based on depth  */
	glEnable(GL_DEPTH_TEST); 
}

static void onDisplay() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glPointSize(5);

	/******* uniform variable initialization*********/
	Matrix4f World;
	World.m[0][0] = cosf(rotation); World.m[0][1] = 0.0f; World.m[0][2] = sinf(rotation); World.m[0][3] = 0.0f;
	World.m[1][0] = 0.0f; World.m[1][1] = 1.0f;  World.m[1][2] = 0.0f; World.m[1][3] = 0.0f;
	World.m[2][0] = -sinf(rotation);        World.m[2][1] = 0.0f;         World.m[2][2] = cosf(rotation); World.m[2][3] = 0.0f;
	World.m[3][0] = 0.0f;        World.m[3][1] = 0.0f;         World.m[3][2] = 0.0f; World.m[3][3] = 1.0f;
	glUniformMatrix4fv(gWorldLocation, 1, GL_TRUE, &World.m[0][0]);
	glUniform3f(gllightColor,1.0f,1.0f,1.0f);
	//glUniform3f(globjectColor,1.0f,1.0f,1.0f);
	glUniform3f(gllightPos,1.0f,5.0f,1.0f);
	glUniform3f(glviewPos,0.0f,0.0f,0.0f);
	glUniform1i(gSampler, 0);
	glUniform1i(gSampler1, 1);
	glUniform1f(gIsovalue, isovalue);
	glUniform1f(gfieldvmin,field->vmin);
	glUniform1f(gfieldvmax,field->vmax);

	glUniform1f(ga,1);
	glUniform1f(gb,1);
	glUniform1f(gc,1);
	glUniform1f(gd,0);
	glUniform1i(side,setside);
	glUniform1i(gslice,enable_slice);

	/********* Draw call*****************************/
	glBindVertexArray(VAO);
	Texture_3d->Bind(GL_TEXTURE0);
	Texture_1d->Bind(GL_TEXTURE1);
	glDrawElements(GL_TRIANGLES, np*3, GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);

	/* check for any errors when rendering */
	GLenum errorCode = glGetError();
	if (errorCode == GL_NO_ERROR) {
		/* double-buffering - swap the back and front buffers */
		glutSwapBuffers();
	} else {
		fprintf(stderr, "OpenGL rendering error %d\n", errorCode);
	}
}

/* pre:  glut window has been resized
 */
static void onReshape(int width, int height) {
	glViewport(0, 0, width, height);
	if (!isFullScreen) {
		theWindowWidth = width;
		theWindowHeight = height;
	}
	// update scene based on new aspect ratio....
}

/* pre:  glut window is not doing anything else
   post: scene is updated and re-rendered if necessary */
static void onIdle() {
	static int oldTime = 0;
	if (isAnimating) {
		int currentTime = glutGet((GLenum) (GLUT_ELAPSED_TIME));
		/* Ensures fairly constant framerate */
		if (currentTime - oldTime > ANIMATION_DELAY) {
			// do animation....
			rotation += 0.005;

			oldTime = currentTime;
			/* compute the frame rate */
			computeFPS();
			/* notify window it has to be repainted */
			glutPostRedisplay();
		}
	}
}

/* pre:  mouse is dragged (i.e., moved while button is pressed) within glut window
   post: scene is updated and re-rendered  */
static void onMouseMotion(int x, int y) {
	/* notify window that it has to be re-rendered */
	glutPostRedisplay();
}

/* pre:  mouse button has been pressed while within glut window
   post: scene is updated and re-rendered */
static void onMouseButtonPress(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		// Left button pressed
	}
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP){
		// Left button un pressed
	}
	/* notify window that it has to be re-rendered */
	glutPostRedisplay();
}

/* pre:  key has been pressed
   post: scene is updated and re-rendered */
static void onAlphaNumericKeyPress(unsigned char key, int x, int y) {
	switch (key) {
			/* toggle animation running */
		case 'a':
			isAnimating = !isAnimating;
			break;
			/* reset */
		case 'r':
			rotation = 0;
			break;
			/* quit! */	
			break;
		case 's':
			if (enable_slice==2) enable_slice=1;
			else enable_slice=2;
			break;
		case 'l':
			setside=1;
			break;
		case 'r':
			setside=2;
			break;
		case 'Q':
		case 'q':
		case 27:
			exit(0);
	}

	/* notify window that it has to be re-rendered */
	glutPostRedisplay();
}

/* pre:  arrow or function key has been pressed
   post: scene is updated and re-rendered */
static void onSpecialKeyPress(int key, int x, int y) {
	/* please do not change function of these keys */
	switch (key) {
			/* toggle full screen mode */
		case GLUT_KEY_F1:
			isFullScreen = !isFullScreen;
			if (isFullScreen) {
				theWindowPositionX = glutGet((GLenum) (GLUT_WINDOW_X));
				theWindowPositionY = glutGet((GLenum) (GLUT_WINDOW_Y));
				glutFullScreen();
			} else {
				glutReshapeWindow(theWindowWidth, theWindowHeight);
				glutPositionWindow(theWindowPositionX, theWindowPositionY);
			}
			break;
		case GLUT_KEY_UP:
			if(isovalue+0.1<=field->vmax)
				isovalue+=0.1;
			break;
		case GLUT_KEY_DOWN:
			if(isovalue-0.1>=field->vmin)
				isovalue-=0.1;
			break;
	}

	/* notify window that it has to be re-rendered */
	glutPostRedisplay();
}

/* pre:  glut window has just been iconified or restored 
   post: if window is visible, animate model, otherwise don't bother */
static void onVisible(int state) {
	if (state == GLUT_VISIBLE) {
		/* tell glut to show model again */
		glutIdleFunc(onIdle);
	} else {
		glutIdleFunc(NULL);
	}
}

static void InitializeGlutCallbacks() {
	/* tell glut how to display model */
	glutDisplayFunc(onDisplay);
	/* tell glutwhat to do when it would otherwise be idle */
	glutIdleFunc(onIdle);
	/* tell glut how to respond to changes in window size */
	glutReshapeFunc(onReshape);
	/* tell glut how to handle changes in window visibility */
	glutVisibilityFunc(onVisible);
	/* tell glut how to handle key presses */
	glutKeyboardFunc(onAlphaNumericKeyPress);
	glutSpecialFunc(onSpecialKeyPress);
	/* tell glut how to handle the mouse */
	glutMotionFunc(onMouseMotion);
	glutMouseFunc(onMouseButtonPress);
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);

	/* request initial window size and position on the screen */
	glutInitWindowSize(theWindowWidth, theWindowHeight);
	glutInitWindowPosition(theWindowPositionX, theWindowPositionY);
	/* request full color with double buffering and depth-based rendering */
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
	/* create window whose title is the name of the executable */
	glutCreateWindow(theProgramTitle);

	InitializeGlutCallbacks();

	// Must be done after glut is initialized!
	GLenum res = glewInit();
	if (res != GLEW_OK) {
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		return 1;
	}

	printf("GL version: %s\n", glGetString(GL_VERSION));

	/* initialize model */
	auto start = high_resolution_clock::now(); 
	onInit(argc, argv);
	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<milliseconds>(stop - start); 
	cout << "Total time=" << duration.count() << endl; 

	/* give control over to glut to handle rendering and interaction  */
	glutMainLoop();

	/* program should never get here */

	return 0;
}

