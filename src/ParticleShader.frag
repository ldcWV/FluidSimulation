#version 330 core

in vec3 normal;
out vec4 FragColor;

void main() {
   vec3 color = vec3(0.5, 0.5, 1.0);
   float dotProd = max(0.0, normal.z);
   FragColor = vec4(dotProd * color, 1.0);
}
