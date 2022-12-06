#version 330 core

in vec3 position;
in vec3 normal;
out vec4 FragColor;

uniform vec3 lightPos;
uniform vec3 cameraPos;

void main() {
   vec3 lightDir = normalize(lightPos - position);
   vec3 viewDir = normalize(cameraPos - position);
   vec3 nnormal = normalize(normal);

   vec3 albedo = vec3(0.5, 0.5, 1.0);
   
   vec3 ambient = 0.1 * albedo;
   vec3 diffuse = 0.5 * albedo * max(0.0, dot(nnormal, lightDir));
   vec3 specular = 0.6 * albedo * pow(max(0.0, dot(reflect(-lightDir, nnormal), viewDir)), 2);

   FragColor = vec4(ambient + diffuse + specular, 1.0);
}
