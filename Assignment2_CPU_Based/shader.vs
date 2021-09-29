#version 130
in vec3 Position;
in vec3 aNormal;
in float infieldval;
in vec3 incolor;
out vec3 Normal;
out vec3 FragPos;
out float fieldval;
out vec3 color;
uniform mat4 gWorld;

void main()
{
    gl_Position = gWorld* vec4(0.6*Position, 1.0);
    FragPos = Position;
    Normal = aNormal;
    fieldval=infieldval;
    color = incolor;
}
