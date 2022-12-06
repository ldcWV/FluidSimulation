#version 330 core

layout (location=0) in vec3 aPos;
layout (location=1) in vec3 translation;

uniform float scale;
uniform mat4 view;
uniform mat4 projection;

out vec3 normal;
out vec3 position;

void main() {
   vec3 newPos = scale * aPos + translation;
   gl_Position = projection * view * vec4(newPos.x, newPos.y, newPos.z, 1.0);
   normal = normalize(aPos);
   position = newPos;
}
