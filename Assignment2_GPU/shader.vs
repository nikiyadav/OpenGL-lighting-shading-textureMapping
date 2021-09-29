#version 130

in vec3 Position;
in vec3 aNormal;
in vec3 texcoor;

out vec3 Normal;
out vec3 FragPos;
out vec3 Texcoor;

uniform mat4 gWorld;

void main()
{
    gl_Position = gWorld* vec4(0.6*Position, 1.0);
    FragPos = Position;
    Normal = aNormal;
    Texcoor=texcoor;
}
