#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b))

typedef struct scalarField {
    float* data;
    float vmin, vmax;
    int dim[3];
    float origin[3];
    float step[3];
}ScalarField;


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
/*
int main(int argc, char ** argv) {
    // testing code
    int i, j;
    ScalarField *field = loadField("2oar.pot");

    printf("Scalar Field\n");
    printf("Dimensions: %d %d %d\n", field->dim[0], field->dim[1], field->dim[2]);
    printf("Origin: %f %f %f\n", field->origin[0], field->origin[1], field->origin[2]);
    printf("Step: %f %f %f\n", field->step[0], field->step[1], field->step[2]);
    printf("Min: %f Max: %f\n", field->vmin, field->vmax);

    int indices[3];
    indices[0] = indices[1] = indices[2] = 0;
    float val = getGridValue(field, indices);
    printf("Value ([%d, %d, %d]): %f\n", indices[0], indices[1], indices[2], val);
    indices[0] = field->dim[0]-1;
    indices[1] = field->dim[1]-1;
    indices[2] = field->dim[2]-1;
    val = getGridValue(field, indices);
    printf("Value ([%d, %d, %d]): %f\n", indices[0], indices[1], indices[2], val);

    freeScalarField(field);
    return 0;
}

*/