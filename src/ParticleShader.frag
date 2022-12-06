#version 330 core

in vec3 normal;
out vec4 FragColor;

uniform vec3 lightDir;

void main() {
   vec3 color = vec3(0.5, 0.5, 1.0);
   float dotProd = max(0.0, dot(normalize(normal), lightDir));
   FragColor = vec4(dotProd * color, 1.0);
}
