#version 330

layout (triangles) in;
layout (triangle_strip, max_vertices=3) out;

in vec3 GPos[];
in vec3 Normal[];
in vec3 Texcoor[];

out vec3 FragPos;
out vec3 FragNormal;
out vec3 FragTexcoor;

uniform float ga, gb, gc, gd;
uniform int side;
uniform int gslice;

int decider(vec4 A, vec4 B, vec4 C) {
    float x1,y1,z1,x2,y2,z2,x3,y3,z3;
    x1=A.x; y1=A.y; z1=A.z;
    x2=B.x; y2=B.y; z2=B.z;
    x3=C.x; y3=C.y; z3=C.z;

    float eval1=ga*x1+gb*y1+gc*z1+gd;
    float eval2=ga*x2+gb*y2+gc*z2+gd;
    float eval3=ga*x3+gb*y3+gc*z3+gd;

    if ( eval1<=0 && eval2<=0 && eval3<=0 ) 
        return 1;
    else if ( eval1>=0 && eval2>=0 && eval3>=0 ) 
        return 2;
    return 3;
}

void main() {

    int l = decider(gl_in[0].gl_Position, gl_in[1].gl_Position, gl_in[2].gl_Position);
    if (gslice==1) {
        if(l == 1 && side == 1) {
            for(int i = 0; i < 3; i++) {
                FragPos = GPos[i];
                FragNormal = Normal[i];
                FragTexcoor = Texcoor[i];
                gl_Position = gl_in[i].gl_Position;
                EmitVertex();
            }
        }
        else if (l==2 && side==2) {
            for(int i = 0; i < 3; i++) {
                FragPos = GPos[i];
                FragNormal = Normal[i];
                FragTexcoor = Texcoor[i];
                gl_Position = gl_in[i].gl_Position;
                EmitVertex();
            }
        }
    }
    else {
        for(int i = 0; i < 3; i++) {
                FragPos = GPos[i];
                FragNormal = Normal[i];
                FragTexcoor = Texcoor[i];
                gl_Position = gl_in[i].gl_Position;
                EmitVertex();
            }
    }
    EndPrimitive();
}
