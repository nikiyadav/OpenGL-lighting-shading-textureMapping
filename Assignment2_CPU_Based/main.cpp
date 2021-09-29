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

GLuint VBO, VAO, nbo;
unsigned int iso_vao,iso_vbo,iso_cbo;
unsigned int slice_vao,slice_vbo,slice_cbo,slice_ibo;
unsigned int grid_vao,grid_vbo,grid_cbo;
unsigned int left_slice_vao,left_slice_vbo,left_slice_cbo,left_slice_nbo;
unsigned int right_slice_vao,right_slice_vbo,right_slice_cbo,right_slice_nbo;

unsigned int left_slice_ver_count, right_slice_ver_count;
int np,n,nv;

float iso_value=0.5;
bool display_isocontour=false;
bool display_left_slice=false, display_right_slice=false;
bool slice=false;

int iso_point_count=0;
int slice_plane_ver_count=0;
int plane_indices_size=0;
GLint gWorldLocation,gllightColor,globjectColor,gllightPos,glviewPos;
//
Vector3f *vertices,*vertices_unnorm,*normals,*scalarColor;
unsigned int *indices;
float *scalarValues;
//

/* Constants */
const int ANIMATION_DELAY = 20; /* milliseconds between rendering */
const char* pVSFileName = "shader.vs";
const char* pFSFileName = "shader.fs";
/********************************************************************
  Utility functions


/****************************************************************/
static float normaliseVertexAttribute(float x, float min, float max) {
	if (max-min)
		return (2*((x-min)/(max-min))-1.0);
	else
		return 0.0f;
} 

static float normaliseScalarField(float s, float smin, float smax ) {
	return (s-smin)/(smax-smin);
}
/*****************************************************************/
/*********************************************************************
POT Reader */

typedef struct scalarField {
    float* data;
    float vmin, vmax;
    int dim[3];
    float origin[3];
    float step[3];
}ScalarField;

ScalarField *field;

/* This function will return the scalar value corresponding to the given grid point. Please ensure that the indices are within the range of grid size (i.e. 0 <= indices[k] < field->dim[k]) */
float getGridValue(ScalarField *field, int *indices){
    int x = indices[0];
    int y = indices[1];
    int z = indices[2];
    int index = z*field->dim[0]*field->dim[1] + y*field->dim[0] + x;
    return field->data[index];
}

/* This function will return the scalar value corresponding to the grid point closest (using floor) to the given point */
float getValue(ScalarField *field, float *point){
    int indices[3];
    if(point[0] < field->origin[0]) point[0] = field->origin[0];
    if(point[1] < field->origin[1]) point[1] = field->origin[1];
    if(point[2] < field->origin[2]) point[2] = field->origin[2];
    indices[0] = min((int)((point[0] - field->origin[0])/ field->step[0]), field->dim[0]-1);
    indices[1] = min((int)((point[1] - field->origin[1])/ field->step[1]), field->dim[1]-1);
    indices[2] = min((int)((point[2] - field->origin[2])/ field->step[2]), field->dim[2]-1);
    return getGridValue(field, indices);
}

float getValueTriLinear(ScalarField *field, float *point){
   // Implement your trilinear interpolation code here
	int indices[3];
    if(point[0] < field->origin[0]) point[0] = field->origin[0];
    if(point[1] < field->origin[1]) point[1] = field->origin[1];
    if(point[2] < field->origin[2]) point[2] = field->origin[2];
    indices[0] = min((int)((point[0] - field->origin[0])/ field->step[0]), field->dim[0]-1);
    indices[1] = min((int)((point[1] - field->origin[1])/ field->step[1]), field->dim[1]-1);
    indices[2] = min((int)((point[2] - field->origin[2])/ field->step[2]), field->dim[2]-1);

    int x=indices[0],y=indices[1],z=indices[2];
    int c000[3]={x,y,z};
    int c100[3]={x+1,y,z};
    int c001[3]={x,y+1,z};
    int c101[3]={x+1,y+1,z};
    int c010[3]={x,y,z+1};
    int c110[3]={x+1,y,z+1};
    int c011[3]={x,y+1,z+1};
    int c111[3]={x+1,y+1,z+1};

    float v1,v2,v3,v4,v5,v6,v7,v8;
    v1 = getGridValue(field,c000);
    v2 = getGridValue(field,c100);
    v3 = getGridValue(field,c001);
    v4 = getGridValue(field,c101);
    v5 = getGridValue(field,c010);
    v6 = getGridValue(field,c110);
    v7 = getGridValue(field,c011);
    v8 = getGridValue(field,c111);

    float xx,yy,zz;
    xx=indices[0]*field->step[0]+field->origin[0];
    yy=indices[1]*field->step[1]+field->origin[1];
    zz=indices[2]*field->step[2]+field->origin[2];

    float p = v1+((v2-v1)/field->step[0])*(point[0]-xx);
    float q = v3+((v4-v3)/field->step[0])*(point[0]-xx);

    float r = p+((q-p)/field->step[1])*(point[1]-yy);

	p = v5+((v6-v5)/field->step[0])*(point[0]-xx);
    q = v7+((v8-v7)/field->step[0])*(point[0]-xx);

    float s = p+((q-p)/field->step[1])*(point[1]-yy);

   	float t = r+((s-r)/field->step[2])*(point[2]-zz);
 
   	return t;
}

ScalarField* loadField(const char* filename){
    ScalarField* field = (ScalarField*) malloc(sizeof(ScalarField));

    FILE* fp = fopen(filename, "r");
    char line[100];
    char temp[40];
    /* Read the size of the grid */
    fgets(line, 100, fp);
    sscanf(line, "%d %d %d", &field->dim[0], &field->dim[1], &field->dim[2]);
    /* Read the origin of the 3D scalar field */
    fgets(line, 100, fp);
    sscanf(line, "%s %f %f %f", temp, &field->origin[0], &field->origin[1], &field->origin[2]);
    /* Read the step size of the grid */
    fgets(line, 100, fp);
    sscanf(line, "%s %f %f %f", temp, &field->step[0], &field->step[1], &field->step[2]);

    /* allocate required data */
    int totalPts = field->dim[0]*field->dim[1]*field->dim[2];
    field->data = (float *) malloc(totalPts * sizeof(float));
    float* tempData = (float *) malloc(totalPts * sizeof(float));
    /* Read the data from file*/
    int i;
    for( i=0;i<totalPts;i++){
        float val;
        fscanf(fp, "%f", &val);
        if (val > 0){
                val = log(1+val);
        } else {
                val = -log(1-val);
        }
        if(i==0){
            field->vmin = field->vmax = val;
        }
        if (val < field->vmin){
            field->vmin = val;
        } else if (val > field->vmax){
            field->vmax = val;
        }
        tempData[i] = val;
    }
    fclose(fp);
    /* Flip order */
    int index = 0,z,y,x;
    for (z=0; z < field->dim[2]; z++){
        for ( y=0; y < field->dim[1]; y++){
            for (x=0; x < field->dim[0]; x++){
                int i = x*field->dim[1]*field->dim[2] + y*field->dim[2] + z;
                field->data[index] = tempData[i];
                index++;
            }
        }
    }
    free(tempData);
    return field;
}

int freeScalarField(ScalarField *field){
    if( field == NULL )
        return 0;
    free(field->data);
    free(field);
    return 1;
}


/*********************************************************************/

/*Slice*/

float t;
static bool HelperToFindIntersectionPoint ( Vector3f X, Vector3f Y , float a, float b, float c, float d) {
	float x1,y1,z1,x2,y2,z2;
	x1 = X.x; y1 = X.y; z1 = X.z;
	x2 = Y.x; y2 = Y.y; z2 = Y.z;

	float deno = a*(x2-x1)+b*(y2-y1)+c*(z2-z1);

	if ( deno ) {
		float num = -1*(a*x1+b*y1+c*z1+d);
		t=num/deno;
		return true;
	}
	else return false;
}

static float EuclideanDistance ( Vector3f X, Vector3f Y ) {
	return sqrt(pow(X.x-Y.x,2)+pow(X.y-Y.y,2)+pow(X.z-Y.z,2));
}

static bool ValidVertex( Vector3f X, Vector3f Y, Vector3f P) {
	float distXY = EuclideanDistance(X,Y);
	float distXP = EuclideanDistance(X,P);
	float distYP = EuclideanDistance(Y,P);

	if ( (distXP > distXY) || (distYP > distXY) ) 
		return false;
	return true;
}

static int decider( Vector3f A,Vector3f B, Vector3f C, float a, float b, float c, float d ) {
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	x1=A.x; y1=A.y; z1=A.z;
	x2=B.x; y2=B.y; z2=B.z;
	x3=C.x; y3=C.y; z3=C.z;

	float eval1=a*x1+b*y1+c*z1+d;
	float eval2=a*x2+b*y2+c*z2+d;
	float eval3=a*x3+b*y3+c*z3+d;

	if ( eval1<=0 && eval2<=0 && eval3<=0 ) 
		return 1;
	else if ( eval1>=0 && eval2>=0 && eval3>=0 ) 
		return 2;
	else 
		return 3;
}

static int deciderForVertex(Vector3f A, float a, float b, float c, float d) {
	float x1,y1,z1;
	x1=A.x; y1=A.y; z1=A.z;
	float eval=a*x1+b*y1+c*z1+d;
	if (eval<=0)
		return 1;
	else 
		return 0;
}

/*Slice model*/
static void sliceModel(float ap, float bp, float cp, float dp) {
	vector<Vector3f> L,R,LC,RC,LN,RN;

	unsigned int a,b,c;
	Vector3f A,B,C;
	for ( int i=0; i<np*3; i+=3 ) {
		a = indices[i]; b = indices[i+1]; c = indices[i+2]; // indices of triangle
		A = vertices_unnorm[a]; B = vertices_unnorm[b]; C = vertices_unnorm[c]; // vertices of triangle
		switch ( decider(A,B,C,ap,bp,cp,dp)) {
			case 1:
				L.push_back(A); LN.push_back(normals[a]);LC.push_back(scalarColor[a]);
				L.push_back(B); LN.push_back(normals[b]);LC.push_back(scalarColor[b]);
				L.push_back(C); LN.push_back(normals[c]);LC.push_back(scalarColor[c]);
				break;
			case 2:
				R.push_back(A); RN.push_back(normals[a]);RC.push_back(scalarColor[a]);
				R.push_back(B); RN.push_back(normals[b]);RC.push_back(scalarColor[b]);
				R.push_back(C); RN.push_back(normals[c]);RC.push_back(scalarColor[c]);
				break;
			/*case 3:
				float x,y,z;
				bool p=false,q=false,s=false;
				Vector3f P,Q,S;
				if (HelperToFindIntersectionPoint(A,B,ap,bp,cp,dp)) {
					x = A.x + t*(B.x-A.x);
					y = A.y + t*(B.y-A.y);
					z = A.z + t*(B.z-A.z);
					P = Vector3f(x,y,z);
					if (ValidVertex(A,B,P)){
						p=true;
					}
				}
				if (HelperToFindIntersectionPoint(B,C,ap,bp,cp,dp)) {
					x = B.x + t*(C.x-B.x);
					y = B.y + t*(C.y-B.y);
					z = B.z + t*(C.z-B.z);
					Q = Vector3f(x,y,z);
					if (ValidVertex(B,C,Q)){
						q=true;
					}
				}
				if (HelperToFindIntersectionPoint(C,A,ap,bp,cp,dp)) {
					x = C.x + t*(A.x-C.x);
					y = C.y + t*(A.y-C.y);
					z = C.z + t*(A.z-C.z);
					S = Vector3f(x,y,z);
					if (ValidVertex(C,A,S)){
						s=true;
					}
				}

				if (p&&q) {

					if ( deciderForVertex(A,ap,bp,cp,dp)) {
						//A.x-=0.5f; P.x-=0.5f; Q.x-=0.5f; C.x-=0.5f;
						L.push_back(A); L.push_back(P); L.push_back(C);
						L.push_back(P); L.push_back(Q); L.push_back(C);
						//P.x+=1.0f; Q.x+=1.0f; B.x+=0.5f;
						R.push_back(P); R.push_back(Q); R.push_back(B);
					}
					else {
						//A.x+=0.5f; P.x+=0.5f; Q.x+=0.5f; C.x+=0.5f;
						R.push_back(A); R.push_back(P); R.push_back(C);
						R.push_back(P); R.push_back(Q); R.push_back(C);
						//P.x-=1.0f; Q.x-=1.0f; B.x-=0.5f;
						L.push_back(P); L.push_back(Q); L.push_back(B);
					}
				}
				else if (q&&s) {
					if ( deciderForVertex(A,ap,bp,cp,dp) ) {
						//A.x-=0.5f; S.x-=0.5f; Q.x-=0.5f; B.x-=0.5f;
						L.push_back(A); L.push_back(S); L.push_back(B);
						L.push_back(S); L.push_back(Q); L.push_back(B);
						//S.x+=1.0f; Q.x+=1.0f; C.x+=0.5f;
						R.push_back(S); R.push_back(Q); R.push_back(C);
					}
					else {
						//A.x+=0.5f; S.x+=0.5f; Q.x+=0.5f; B.x+=0.5f;
						R.push_back(A); R.push_back(S); R.push_back(B);
						R.push_back(S); R.push_back(Q); R.push_back(B);
						//S.x-=1.0f; Q.x-=1.0f; C.x-=0.5f;
						L.push_back(S); L.push_back(Q); L.push_back(C);
					}

				}
				else if (p&&s) {
					if ( deciderForVertex(B,ap,bp,cp,dp) ) {
						//B.x-=0.5f; S.x-=0.5f; P.x-=0.5f; C.x-=0.5f;
						L.push_back(B); L.push_back(P); L.push_back(C);
						L.push_back(P); L.push_back(S); L.push_back(C);
						//S.x+=1.0f; P.x+=1.0f; A.x+=0.5f;
						R.push_back(S); R.push_back(P); R.push_back(A);
					}
					else {
						//B.x+=0.5f; S.x+=0.5f; P.x+=0.5f; C.x+=0.5f;
						R.push_back(B); R.push_back(P); R.push_back(S);
						R.push_back(B); R.push_back(C); R.push_back(S);
						//S.x-=1.0f; P.x-=1.0f; A.x-=0.5f;
						L.push_back(S); L.push_back(P); L.push_back(A);
					}
				}*/

		}
	}

	left_slice_ver_count = L.size();
	right_slice_ver_count = R.size();

	Vector3f *left_color_buffer = (Vector3f *)malloc(sizeof(Vector3f)*L.size());
	Vector3f *left_normals = (Vector3f *)malloc(sizeof(Vector3f)*L.size());
	for ( int i=0;i<L.size();i++ ) {
		left_color_buffer[i]=LC[i];
		left_normals[i]=LN[i];
	}

	Vector3f *right_color_buffer = (Vector3f *)malloc(sizeof(Vector3f)*R.size());
	Vector3f *right_normals = (Vector3f *)malloc(sizeof(Vector3f)*R.size());
	for ( int i=0;i<R.size();i++ ) {
		right_color_buffer[i]=RC[i];
		right_normals[i]=RN[i];
	}

	for ( int i=0;i<L.size();i++ ) {
		L[i].x=normaliseVertexAttribute(L[i].x,field->origin[0],field->origin[0]+(field->dim[0]-1)*field->step[0]);
		L[i].y=normaliseVertexAttribute(L[i].y,field->origin[1],field->origin[1]+(field->dim[1]-1)*field->step[1]);
		L[i].z=normaliseVertexAttribute(L[i].z,field->origin[2],field->origin[2]+(field->dim[2]-1)*field->step[2]);
	}

	for ( int i=0;i<R.size();i++ ) {
		R[i].x=normaliseVertexAttribute(R[i].x,field->origin[0],field->origin[0]+(field->dim[0]-1)*field->step[0]);
		R[i].y=normaliseVertexAttribute(R[i].y,field->origin[1],field->origin[1]+(field->dim[1]-1)*field->step[1]);
		R[i].z=normaliseVertexAttribute(R[i].z,field->origin[2],field->origin[2]+(field->dim[2]-1)*field->step[2]);
	}

	glGenVertexArrays(1, &left_slice_vao);
	glGenVertexArrays(1,&right_slice_vao);

	// Left
	glBindVertexArray(left_slice_vao);

	glGenBuffers(1,&left_slice_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, left_slice_vbo);
	glBufferData(GL_ARRAY_BUFFER, L.size()*sizeof(Vector3f), &L[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE, 0, 0);

	glGenBuffers(1,&left_slice_nbo);
	glBindBuffer(GL_ARRAY_BUFFER, left_slice_nbo);
	glBufferData(GL_ARRAY_BUFFER, L.size()*sizeof(Vector3f), left_normals, GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE, 0, 0);

	glGenBuffers(1,&left_slice_cbo);
	glBindBuffer(GL_ARRAY_BUFFER, left_slice_cbo);
	glBufferData(GL_ARRAY_BUFFER, L.size()*sizeof(Vector3f), left_color_buffer, GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE, 0, 0);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER,0);

	//Right
	glBindVertexArray(right_slice_vao);

	glGenBuffers(1,&right_slice_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, right_slice_vbo);
	glBufferData(GL_ARRAY_BUFFER, R.size()*sizeof(Vector3f), &R[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE, 0, 0);

	glGenBuffers(1,&right_slice_nbo);
	glBindBuffer(GL_ARRAY_BUFFER, right_slice_nbo);
	glBufferData(GL_ARRAY_BUFFER, R.size()*sizeof(Vector3f), right_normals, GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE, 0, 0);

	glGenBuffers(1,&right_slice_cbo);
	glBindBuffer(GL_ARRAY_BUFFER, right_slice_cbo);
	glBufferData(GL_ARRAY_BUFFER, R.size()*sizeof(Vector3f), right_color_buffer, GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE, 0, 0);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER,0);
}

/*slice grid*/
// plane equation ax+by+cz+d=0
static void sliceGrid(float a, float b, float c, float d) {
	vector<Vector3f> IV; 
	std::vector<unsigned int> PlaneIndices;

	int i,j,k;
	for (i=0;i<field->dim[0];i++) {
		for (j=0;j<field->dim[1];j++) {
			for (k=0;k<field->dim[2];k++) {
				float x1,y1,z1,x2,y2,z2;
				//bottom left voxel point
				x1 = field->origin[0]+i*field->step[0];
				y1 = field->origin[1]+j*field->step[1];
				z1 = field->origin[2]+k*field->step[2];

				//top right voxel point
				x2 = x1+field->step[0];
				y2 = y1+field->step[1];
				z2 = z1+field->step[2];

				bool intersect[12];
				for ( int l=0;l<12;l++ ) intersect[l]=false;
				Vector3f intersectionPoints[12];
				for ( int l=0;l<12;l++ ) intersectionPoints[l]=Vector3f(0,0,0);
				float x, y, z;

				if (HelperToFindIntersectionPoint( Vector3f(x1,y1,z1), Vector3f(x2,y1,z1), a, b, c, d )) {
					x = x1+(x2-x1)*t;
					y = y1;
					z = z1;
					intersectionPoints[0] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x1,y1,z1), Vector3f(x2,y1,z1), intersectionPoints[0]) ) 
						intersect[0]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x1,y2,z1), Vector3f(x2,y2,z1), a, b, c, d )) {
					x = x1+(x2-x1)*t;
					y = y2;
					z = z1;
					intersectionPoints[1] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x1,y2,z1), Vector3f(x2,y2,z1), intersectionPoints[1]) ) 
						intersect[1]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x1,y1,z2), Vector3f(x2,y1,z2), a, b, c, d )) {
					x = x1+(x2-x1)*t;
					y = y1;
					z = z2;
					intersectionPoints[2] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x1,y1,z2), Vector3f(x2,y1,z2), intersectionPoints[2]) ) 
						intersect[2]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x1,y2,z2), Vector3f(x2,y2,z2), a, b, c, d )) {
					x = x1+(x2-x1)*t;
					y = y2;
					z = z2;
					intersectionPoints[3] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x1,y2,z2), Vector3f(x2,y2,z2), intersectionPoints[3] ) ) 
						intersect[3]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x1,y1,z1), Vector3f(x1,y1,z2), a, b, c, d )) {
					x = x1;
					y = y1;
					z = z1+(z2-z1)*t;
					intersectionPoints[4] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x1,y1,z1), Vector3f(x1,y1,z2), intersectionPoints[4] ) ) 
						intersect[4]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x1,y2,z1), Vector3f(x1,y2,z2), a, b, c, d )) {
					x = x1;
					y = y2;
					z = z1+(z2-z1)*t;
					intersectionPoints[5] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x1,y2,z1), Vector3f(x1,y2,z2), intersectionPoints[5] ))
						intersect[5]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x2,y1,z1), Vector3f(x2,y1,z2), a, b, c, d )) {
					x = x2;
					y = y1;
					z = z1+(z2-z1)*t;
					intersectionPoints[6] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x2,y1,z1), Vector3f(x2,y1,z2), intersectionPoints[6] ) ) 
						intersect[6]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x2,y2,z1), Vector3f(x2,y2,z2), a, b, c, d )) {
					x = x2;
					y = y2;
					z = z1+(z2-z1)*t;
					intersectionPoints[7] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x2,y2,z1), Vector3f(x2,y2,z2), intersectionPoints[7] )) 
						intersect[7]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x1,y1,z1), Vector3f(x1,y2,z1), a, b, c, d )) {
					x = x1;
					y = y1+(y2-y1)*t;
					z = z1;
					intersectionPoints[8]= Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x1,y1,z1), Vector3f(x1,y2,z1), intersectionPoints[8] ))  
						intersect[8]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x1,y1,z2), Vector3f(x1,y2,z2), a, b, c, d )) {
					x = x1;
					y = y1+(y2-y1)*t;
					z = z2;
					intersectionPoints[9] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x1,y1,z2), Vector3f(x1,y2,z2), intersectionPoints[9] )) 
						intersect[9]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x2,y1,z1), Vector3f(x2,y2,z1), a, b, c, d )) {
					x = x2;
					y = y1+(y2-y1)*t;
					z = z1;
					intersectionPoints[10] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x2,y1,z1), Vector3f(x2,y2,z1), intersectionPoints[10] )) 
						intersect[10]=true;
				}
				if (HelperToFindIntersectionPoint( Vector3f(x2,y1,z2), Vector3f(x2,y2,z2), a, b, c, d )) {
					x = x2;
					y = y1+(y2-y1)*t;
					z = z2;
					intersectionPoints[11] = Vector3f(x,y,z);
					if ( ValidVertex( Vector3f(x2,y1,z2), Vector3f(x2,y2,z2), intersectionPoints[11] )) 
						intersect[11]=true;
				}

				int count=0;
				for (int l=0;l<12;l++) {
					if (intersect[l])
						count++; 
				}

				unsigned int size=0;
				switch(count) {
					case 0:
						// no intersection of voxel with plane
						break;
					case 1:
						// plane intersect voxel at a single point
						// Add this point to vbo 3 times for a triangle
						size=IV.size();
						for ( int l=0;l<12;l++ ) {
							if ( intersect[l]) {
								IV.push_back(intersectionPoints[l]);
							}
						}
						PlaneIndices.push_back(size);
						PlaneIndices.push_back(size);
						PlaneIndices.push_back(size);
						break;
					case 3:
						// Add these 3 points to vbo
						size=IV.size();
						for ( int l=0;l<12;l++ ) {
							if (intersect[l]) {
								IV.push_back(intersectionPoints[l]);
							}
						}

						PlaneIndices.push_back(size); PlaneIndices.push_back(size+1); PlaneIndices.push_back(size+2);
						break;
					case 4:
						size=IV.size();
						for ( int l=0;l<12;l++ ) {
							if (intersect[l]) {
								IV.push_back(intersectionPoints[l]);
							}
						}
						// Add all possible 4 traingles

						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+3);
						PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);
						break;
					case 5:
						size=IV.size();
						for ( int l=0;l<12;l++ ) {
							if (intersect[l]) {
								IV.push_back(intersectionPoints[l]);
							}
						}

						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+3);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);
						PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+4);
						break;
					case 6:
						size=IV.size();
						for ( int l=0;l<12;l++ ) {
							if (intersect[l]) {
								IV.push_back(intersectionPoints[l]);
							}
						}
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+3);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+5);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+4);PlaneIndices.push_back(size+5);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+5);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+5);
						PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);
						PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+5);
						PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+5);
						PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+4);PlaneIndices.push_back(size+5);
						PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+4);PlaneIndices.push_back(size+1);
						PlaneIndices.push_back(size+4);PlaneIndices.push_back(size+5);PlaneIndices.push_back(size+1);
						PlaneIndices.push_back(size+4);PlaneIndices.push_back(size+5);PlaneIndices.push_back(size+2);
						PlaneIndices.push_back(size+2);PlaneIndices.push_back(size+3);PlaneIndices.push_back(size+4);
						PlaneIndices.push_back(size+0);PlaneIndices.push_back(size+1);PlaneIndices.push_back(size+2); // change this
				}
			}
		}
	}

	slice_plane_ver_count=IV.size();
	plane_indices_size=PlaneIndices.size();
	//cout << "slice_plane_ver_count=" << slice_plane_ver_count << "\n";

	/************************finding scalar values on plane vertices***********/
	float *scalarValuesPlane=(float *)malloc(sizeof(float)*IV.size());

	for (i=0;i<IV.size();i++) {
		float vertex[3];
		vertex[0]=IV[i].x; vertex[1]=IV[i].y; vertex[2]=IV[i].z;		
		float val = getValueTriLinear(field,vertex);
		val = normaliseScalarField(val,field->vmin, field->vmax);
		scalarValuesPlane[i]=val;
	}

	/*********************color values on plane vertices***************************/
	Vector3f *slice_plane_color=(Vector3f *)malloc(sizeof(Vector3f)*IV.size());

	Vector3f Red=Vector3f(1.0f,0.0f,0.0f);
	Vector3f Green=Vector3f(0.0f,0.0f,1.0f);
	for (int i=0;i<IV.size();i++) {
		slice_plane_color[i] = Red + (Green-Red)*scalarValuesPlane[i];
	}

	for ( int i=0;i<IV.size();i++ ) {
		IV[i].x=normaliseVertexAttribute(IV[i].x,field->origin[0],field->origin[0]+(field->dim[0]-1)*field->step[0]);
		IV[i].y=normaliseVertexAttribute(IV[i].y,field->origin[1],field->origin[1]+(field->dim[1]-1)*field->step[1]);
		IV[i].z=normaliseVertexAttribute(IV[i].z,field->origin[2],field->origin[2]+(field->dim[2]-1)*field->step[2]);
	} 

	glBindVertexArray(slice_vao);

	glGenBuffers(1,&slice_vbo);
	glBindBuffer(GL_ARRAY_BUFFER,slice_vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vector3f)*IV.size(),&IV[0],GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);

	glGenBuffers(1,&slice_cbo);
	glBindBuffer(GL_ARRAY_BUFFER,slice_cbo);
	glBufferData(GL_ARRAY_BUFFER,sizeof(Vector3f)*IV.size(),slice_plane_color,GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE,0,0);

	glGenBuffers(1,&slice_ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, slice_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*PlaneIndices.size(), &PlaneIndices[0], GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

}
/********************************************************************/

/*******************************************************************/
//Display Grid

static void displayGridCube() {
	float x1,y1,z1,x2,y2,z2;

	//bottom left voxel point
	x1 = field->origin[0];
	y1 = field->origin[1];
	z1 = field->origin[2];

	//top right voxel point
	x2 = x1+(field->dim[0]-1)*field->step[0];
	y2 = y1+(field->dim[1]-1)*field->step[1];
	z2 = z1+(field->dim[2]-1)*field->step[2];

	Vector3f GridEdges[24];
	GridEdges[0]=Vector3f(x1,y1,z1); GridEdges[1]=Vector3f(x2,y1,z1); 
	GridEdges[2]=Vector3f(x1,y2,z1); GridEdges[3]=Vector3f(x2,y2,z1); 
	GridEdges[4]=Vector3f(x1,y1,z2); GridEdges[5]=Vector3f(x2,y1,z2); 
	GridEdges[6]=Vector3f(x1,y2,z2); GridEdges[7]=Vector3f(x2,y2,z2); 
	GridEdges[8]=Vector3f(x1,y1,z1); GridEdges[9]=Vector3f(x1,y1,z2); 
	GridEdges[10]=Vector3f(x1,y2,z1); GridEdges[11]=Vector3f(x1,y2,z2); 
	GridEdges[12]=Vector3f(x2,y1,z1); GridEdges[13]=Vector3f(x2,y1,z2); 
	GridEdges[14]=Vector3f(x2,y2,z1); GridEdges[15]=Vector3f(x2,y2,z2); 
	GridEdges[16]=Vector3f(x1,y1,z1); GridEdges[17]=Vector3f(x1,y2,z1); 
	GridEdges[18]=Vector3f(x1,y1,z2); GridEdges[19]=Vector3f(x1,y2,z2); 
	GridEdges[20]=Vector3f(x2,y1,z1); GridEdges[21]=Vector3f(x2,y2,z1); 
	GridEdges[22]=Vector3f(x2,y1,z2); GridEdges[23]=Vector3f(x2,y2,z2); 

	Vector3f EdgeColor[24];

	for ( int i=0;i<24;i++ ) {
		EdgeColor[i]=Vector3f(1.0f,1.0f,1.0f);
	}

	for (int i=0;i<24;i++ ) {
		GridEdges[i].x=normaliseVertexAttribute(GridEdges[i].x,x1,x2);
		GridEdges[i].y=normaliseVertexAttribute(GridEdges[i].y,y1,y2);
		GridEdges[i].z=normaliseVertexAttribute(GridEdges[i].z,z1,z2);
		//cout << GridEdges[i].x << " " << GridEdges[i].y << " " << GridEdges[i].z << "\n";
	} 
	
	glBindVertexArray(grid_vao);

	glGenBuffers(1,&grid_vbo);
	glBindBuffer(GL_ARRAY_BUFFER,grid_vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vector3f)*24,GridEdges,GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);

	glGenBuffers(1,&grid_cbo);
	glBindBuffer(GL_ARRAY_BUFFER,grid_cbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vector3f)*24,EdgeColor,GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE,0,0);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/******************************************************************/

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


static void CreateVertexBuffer() {
	glGenVertexArrays(1, &VAO);
	glGenVertexArrays(1,&iso_vao);
	glGenVertexArrays(1,&slice_vao);
	glGenVertexArrays(1,&grid_vao);
	glBindVertexArray(VAO);

	const char* OffFile = "1grm.off";
	int noEdges,i,j;
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

	vertices=(Vector3f *)malloc(sizeof(Vector3f)*nv);
	vertices_unnorm=(Vector3f *)malloc(sizeof(Vector3f)*nv);
	float x,y,z;
	/* Read the vertices' location*/	
	for(i = 0;i < nv;i++) {
		fscanf(input, "%f %f %f", &x,&y,&z);
		vertices[i].x = x;
		vertices[i].y = y;
		vertices[i].z = z;
		vertices_unnorm[i].x = x;
		vertices_unnorm[i].y = y;
		vertices_unnorm[i].z = z;
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

	scalarValues=(float *)malloc(sizeof(float)*nv);

	for (i=0;i<nv;i++) {
		float vertex[3];
		vertex[0]=vertices[i].x; vertex[1]=vertices[i].y; vertex[2]=vertices[i].z;		
		float val = getValueTriLinear(field,vertex);
		val = normaliseScalarField(val,field->vmin, field->vmax);
		scalarValues[i]=val;
	}

	printf("Scalar Field\n");
    printf("Dimensions: %d %d %d\n", field->dim[0], field->dim[1], field->dim[2]);
    printf("Origin: %f %f %f\n", field->origin[0], field->origin[1], field->origin[2]);
    printf("Step: %f %f %f\n", field->step[0], field->step[1], field->step[2]);
    printf("Min: %f Max: %f\n", field->vmin, field->vmax);

	/*Normalize vertices*/
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
	indices=(unsigned int*)malloc(sizeof(unsigned int)*np*3);
	j=0;

	for( std::vector<unsigned int>::iterator it=faces.begin(); it!=faces.end(); it++ ){
		indices[j]=*it;
		j++;
	}

	/* Calculate normals*/
	normals=(Vector3f *)malloc(sizeof(Vector3f)*nv);

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

	//Color Buffer
	scalarColor=(Vector3f *)malloc(sizeof(Vector3f)*nv);

	Vector3f Red=Vector3f(1.0f,0.0f,0.0f);
	Vector3f Blue=Vector3f(0.0f,0.0f,1.0f);
	for (int i=0;i<nv;i++) {
		scalarColor[i] = Red + (Blue-Red)*scalarValues[i];
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

	unsigned int sbo;
	glGenBuffers(1,&sbo);
	glBindBuffer(GL_ARRAY_BUFFER,sbo);
	glBufferData(GL_ARRAY_BUFFER,sizeof(float)*nv,scalarValues,GL_STATIC_DRAW);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,0,0);

	unsigned int cbo;
	glGenBuffers(1,&cbo);
	glBindBuffer(GL_ARRAY_BUFFER,cbo);
	glBufferData(GL_ARRAY_BUFFER,sizeof(Vector3f)*nv,scalarColor,GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE,0,0);

	unsigned int ibo;
	glGenBuffers(1,&ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*np*sizeof(unsigned int), indices, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

}

static void CreateIsoContour( float iso_value ) {
	std::vector<Vector3f> isoPoints;
	int limit=np*n;
	for ( int i=0;i<limit;i+=3 ) {
		// extract a traingle
		int idx0=indices[i];
		int idx1=indices[i+1];
		int idx2=indices[i+2];
		// find scalar values at each vertex of triangle
		float a,b,c;
		a=scalarValues[idx0];
		b=scalarValues[idx1];
		c=scalarValues[idx2]; 
		// find the coordinates of the traingle points.
		Vector3f A,B,C;
		A=vertices[idx0];
		B=vertices[idx1];
		C=vertices[idx2];

		if ( (a<iso_value && b<iso_value && c<iso_value) || (a>iso_value && b>iso_value && c>iso_value) ) {
			continue;
		}
		else {
			//label points: 0 if =iso_value, 1 if <iso_value, 2 if >iso_value
			int la=0,lb=0,lc=0;
			float d;
			d=a-iso_value;
			if (d<0) la=1;
			else if (d>0) la=2;

			d=b-iso_value;
			if (d<0) lb=1;
			else if (d>0) lb=2;

			d=c-iso_value;
			if (d<0) lc=1;
			else if (d>0) lc=2;

			Vector3f U,V;
			// one iso_point is one of the vertex of triangle
			if (la==0) {
				if ( (lb==1 && lc==2) || (lb==2 && lc==1) ) {
					Vector3f W = B + (C-B)*((iso_value-b)/(c-b));
					U=A; V=W;
					isoPoints.push_back(U);
					isoPoints.push_back(V);
				}
				else {
					U=A; V=A;
					isoPoints.push_back(U);
					isoPoints.push_back(V);
				}
			}
			if (lb==0) {
				if ( (la==1 && lc==2) || (la==2 && lc==1) ) {
					Vector3f W = A + (C-A)*((iso_value-a)/(c-a));
					U=B; V=W;
					isoPoints.push_back(U);
					isoPoints.push_back(V);
				}
				else {
					U=B; V=B;
					isoPoints.push_back(U);
					isoPoints.push_back(V);
				}
			}
			if (lc==0) {
				if ( (la==1 && lb==2) || (la==2 && lb==1) ) {
					Vector3f W = A + (B-A)*((iso_value-a)/(b-a));
					U=C; V=W;
					isoPoints.push_back(U);
					isoPoints.push_back(V);
				}
				else {
					U=C; V=C;
					isoPoints.push_back(U);
					isoPoints.push_back(V);
				}
			}
			// two plus points or two minus points exist 
			else {
				// two plus points
				if ( (la==2 && lb==2) || (lb==2 && lc==2) || (lc==2 && la==2)) {
					if (la==1) {
						U = A + (B-A)*((iso_value-a)/(b-a));
						V = A + (C-A)*((iso_value-a)/(c-a));
					}
					else if (lb==1) {
						U = B + (A-B)*((iso_value-b)/(a-b));
						V = B + (C-B)*((iso_value-b)/(c-b)); 
					}
					else if (lc==1) {
						U = C + (A-C)*((iso_value-c)/(a-c));
						V = C + (B-C)*((iso_value-c)/(b-c));
					}
				}
				// two minus points
				else if ( (la==1 && lb==1) || (lb==1 && lc==1) || (lc==1 && la==1)) {
					if (la==2) {
						U = A + (B-A)*((iso_value-a)/(b-a));
						V = A + (C-A)*((iso_value-a)/(c-a));
					}
					else if (lb==2) {
						U = B + (A-B)*((iso_value-b)/(a-b));
						V = B + (C-B)*((iso_value-b)/(c-b)); 
					}
					else if (lc==2) {
						U = C + (A-C)*((iso_value-c)/(a-c));
						V = C + (B-C)*((iso_value-c)/(b-c));
					}
				}
				isoPoints.push_back(U);
				isoPoints.push_back(V);
			}
		}
	}

	iso_point_count=isoPoints.size();

	Vector3f *iso_line_color=(Vector3f *)malloc(sizeof(Vector3f)*iso_point_count);

	for(int i=0;i<iso_point_count;i++) {
		iso_line_color[i]=Vector3f(0.0f,0.0f,0.0f);
	}

	glBindVertexArray(iso_vao);	
	glGenBuffers(1, &iso_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, iso_vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vector3f)*isoPoints.size(), &isoPoints[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glGenBuffers(1,&iso_cbo);
	glBindBuffer(GL_ARRAY_BUFFER,iso_cbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vector3f)*isoPoints.size(), iso_line_color, GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE,0,0);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
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

	glValidateProgram(ShaderProgram);
	glGetProgramiv(ShaderProgram, GL_VALIDATE_STATUS, &Success);
	if (!Success) {
		glGetProgramInfoLog(ShaderProgram, sizeof (ErrorLog), NULL, ErrorLog);
		fprintf(stderr, "Invalid shader program: '%s'\n", ErrorLog);
		exit(1);
	}

	glUseProgram(ShaderProgram);
	gWorldLocation = glGetUniformLocation(ShaderProgram, "gWorld");
	gllightColor = glGetUniformLocation(ShaderProgram,"lightColor");
	//globjectColor = glGetUniformLocation(ShaderProgram,"objectColor");
	gllightPos = glGetUniformLocation(ShaderProgram,"lightPos");
	glviewPos = glGetUniformLocation(ShaderProgram,"viewPos");
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
	//auto start = high_resolution_clock::now(); 
	CreateVertexBuffer();
	//auto stop = high_resolution_clock::now(); 
	//auto duration = duration_cast<milliseconds>(stop - start); 
	//cout << "CreateVertexBuffer=" << duration.count() << endl; 
	CompileShaders();
/*	sliceGrid(-1,1.12,-0.81,20.3);
	sliceModel(-1,1.12,-0.81,20.3);	
	CreateIsoContour(iso_value);
*/
/*	slice=true;
	display_isocontour=false;
	display_right_slice=false;
	display_left_slice=false;
	auto start = high_resolution_clock::now(); 
	sliceGrid(-1,1.12,-0.81,20.3);
	sliceModel(-1,1.12,-0.81,20.3);	
	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<milliseconds>(stop - start); 
	cout << "slice=" << duration.count() << endl; 
	displayGridCube();*/

/*	auto start = high_resolution_clock::now(); 
	CreateIsoContour(iso_value);
	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<milliseconds>(stop - start); 
	cout << "CreateIsoContour=" << duration.count() << endl; 
*/
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

	/********* Draw call*****************************/
	if (!slice) {
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, np*3, GL_UNSIGNED_INT, nullptr);
		glBindVertexArray(0);

		if (display_isocontour) {
			glBindVertexArray(iso_vao);
			glLineWidth(2.5);
			glDrawArrays(GL_LINES, 0, iso_point_count);
			glBindVertexArray(0);
		}
	}
	else {
		//Display Bounding Box
		glBindVertexArray(grid_vao);
		glLineWidth(2.5);
		glDrawArrays(GL_LINES, 0, 24);
		glBindVertexArray(0);

		//Display Plane
		glBindVertexArray(slice_vao);
		glDrawElements(GL_TRIANGLES, plane_indices_size, GL_UNSIGNED_INT, nullptr);
		glBindVertexArray(0);
		
		//Display Left slice
		if (display_left_slice) {
			glBindVertexArray(left_slice_vao);
			glDrawArrays(GL_TRIANGLES,0,left_slice_ver_count);
			glBindVertexArray(0);
		}
		
		//Display Right slice
		if (display_right_slice) {
			glBindVertexArray(right_slice_vao);
			glDrawArrays(GL_TRIANGLES,0,right_slice_ver_count);
			glBindVertexArray(0);
		}
	}

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
	std::chrono::_V2::system_clock::time_point start, stop;
	std::chrono:: duration<long int, std::ratio<1,1000> > duration;
	switch (key) {
			/* toggle animation running */
		case 'a':
			isAnimating = !isAnimating;
			break;
			/* reset */
		case 'r':
			rotation = 0;
			break;
		case 'i':
			display_isocontour = !display_isocontour;
			iso_value=0.5;
			break;
		case 's':
			slice=!slice;
			break;
			/* quit! */
		case 'z':
			slice=true;
			display_isocontour=false;
			display_right_slice=false;
			display_left_slice=false; 
			start = high_resolution_clock::now(); 
			sliceGrid(-1,1.12,-0.81,20.3);
			sliceModel(-1,1.12,-0.81,20.3);	
			displayGridCube();
			stop = high_resolution_clock::now(); 
			duration = duration_cast<milliseconds>(stop - start); 
			cout << "Time taken in slicing=" << duration.count() << endl; 
			break;
		case 'x':
			slice=true;
			display_isocontour=false;
			display_right_slice=false;
			display_left_slice=false;
			start = high_resolution_clock::now(); 
			sliceGrid(0,0,5,0);
			sliceModel(0,0,5,0);
			displayGridCube();
			stop = high_resolution_clock::now(); 
			duration = duration_cast<milliseconds>(stop - start); 
			cout << "Time taken in slicing=" << duration.count() << endl; 
			break;
		case 'c':
			slice=true;
			display_isocontour=false;
			display_right_slice=false;
			display_left_slice=false;
			start = high_resolution_clock::now(); 
			sliceGrid(-1,1,0,5);
			sliceModel(-1,1,0,5);
			displayGridCube();
			stop = high_resolution_clock::now(); 
			duration = duration_cast<milliseconds>(stop - start); 
			cout << "Time taken in slicing=" << duration.count() << endl; 
			break;
		case 'v':
			slice = true;
			display_isocontour=false;
			display_right_slice=false;
			display_left_slice=false;
			start = high_resolution_clock::now(); 
			sliceGrid(1,2,1,15);
			sliceModel(1,2,1,15);
			displayGridCube();
			stop = high_resolution_clock::now(); 
			duration = duration_cast<milliseconds>(stop - start); 
			cout << "Time taken in slicing=" << duration.count() << endl; 
			break;
		case 'b':
			slice = true;
			display_isocontour=false;
			display_right_slice=false;
			display_left_slice=false;
			start = high_resolution_clock::now(); 
			sliceGrid(1,1,1,0);
			sliceModel(1,1,1,0);
			displayGridCube();
			stop = high_resolution_clock::now(); 
			duration = duration_cast<milliseconds>(stop - start); 
			cout << "Time taken in slicing=" << duration.count() << endl; 
			break;
		case 'n':
			display_left_slice=true;
			display_right_slice=false;
			break;
		case 'm':
			display_right_slice=true;
			display_left_slice=false;
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
	/*std::chrono::_V2::system_clock::time_point start, stop;
	std::chrono:: duration<long int, std::ratio<1,1000> > duration;*/
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
			if (display_isocontour) {
				if ((iso_value+0.01)<=1.0f)
					iso_value+=0.01;
				auto start = high_resolution_clock::now(); 
				CreateIsoContour(iso_value);
				auto stop = high_resolution_clock::now(); 
				auto duration = duration_cast<milliseconds>(stop - start); 
				cout << "Time taken to create IsoContour=" << duration.count() << endl; 
			}
			break;
		case GLUT_KEY_DOWN:
			if (display_isocontour) {
				if ((iso_value-0.01)>0.0f)
					iso_value-=0.01;
				auto start = high_resolution_clock::now(); 
				CreateIsoContour(iso_value);
				auto stop = high_resolution_clock::now(); 
				auto duration = duration_cast<milliseconds>(stop - start); 
				cout << "Time taken to create IsoContour=" << duration.count() << endl; 
			}
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
	//auto start = high_resolution_clock::now(); 
	onInit(argc, argv);
	//auto stop = high_resolution_clock::now(); 
	//auto duration = duration_cast<milliseconds>(stop - start); 
	//cout << "Total Time=" << duration.count() << endl; 


	/* give control over to glut to handle rendering and interaction  */
	glutMainLoop();

	/* program should never get here */

	return 0;
}

