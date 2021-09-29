#version 130

out vec4 FragColor;
in vec3 FragPos;
in vec3 Normal;
in float fieldval;
in vec3 color;
uniform vec3 lightColor;
//uniform vec3 objectColor;
uniform vec3 lightPos;
uniform vec3 viewPos;

void main()
{
	//vec3 objectColor = mix(vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), fieldval);
	vec3 objectColor = color;
	//Ambient 
	float ambientStrength = 0.8;
	vec3 ambient=ambientStrength*lightColor;

	//Diffuse
	vec3 norm = normalize(Normal);
	vec3 lightDir =  normalize(lightPos-FragPos);
	float diff = max(dot(norm,lightDir),0.0);
	vec3 diffuse = diff*lightColor;

	//Specular
	float specularStrength=0.5;
	vec3 viewDir = normalize(viewPos-FragPos);
	vec3 reflectDir = reflect(-lightDir, norm);
	float spec = pow(max(dot(viewDir,reflectDir),0.0),5000);
	vec3 specular = specularStrength * spec * lightColor;

	vec3 result=(ambient+diffuse+specular)*objectColor;
	FragColor=vec4(result,1.0f);
	//FragColor=vec4(Normal,1.0f);
	//FragColor=vec4(objectColor,1.0f);

	//float Value = gl_FragCoord.z;
	//outputColor = mix(vec4(1.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f, 0.0f, 1.0f, 1.0f), Value);
}
