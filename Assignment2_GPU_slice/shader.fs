#version 330

out vec4 FragColor;

in vec3 FragPos;
in vec3 FragNormal;
in vec3 FragTexcoor;

uniform sampler3D gSampler;
uniform sampler1D gSampler1;
uniform float gIsovalue;
uniform float fieldvmin;
uniform float fieldvmax;

uniform vec3 lightColor;
//uniform vec3 objectColor;
uniform vec3 lightPos;
uniform vec3 viewPos;

void main()
{
   //Ambient 
    float ambientStrength = 0.8;
    vec3 ambient=ambientStrength*lightColor;

    //Diffuse
    vec3 norm = normalize(FragNormal);
    vec3 lightDir =  normalize(lightPos-FragPos);
    float diff = max(dot(norm,lightDir),0.0);
    vec3 diffuse = diff*lightColor;

    //Specular
    float specularStrength=0.5;
    vec3 viewDir = normalize(viewPos-FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir,reflectDir),0.0),5000);
    vec3 specular = specularStrength * spec * lightColor;

    vec3 result=(ambient+diffuse+specular);
    FragColor=vec4(result,1.0f);
    //FragColor=vec4(Normal,1.0f);
    //FragColor=vec4(objectColor,1.0f);

    float alpha=0.05;
    vec4 scalar = texture3D(gSampler, FragTexcoor);        
    
    float scalarvalue = (scalar.x-fieldvmin)/(fieldvmax-fieldvmin);
     
    vec4 colormap=texture1D(gSampler1,scalarvalue);
    
    if( scalar.x <= gIsovalue+alpha && scalar.x >= gIsovalue-alpha )
        FragColor =vec4(0.0,0.0,0.0,1.0);
    else
        FragColor =vec4(result,1.0f)*colormap;
    //FragColor=vec4(1.0f,0.0f,0.0f,1.0f);
}
