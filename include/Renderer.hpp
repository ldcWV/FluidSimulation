#pragma once
#include "Scene.hpp"
#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include "Camera.hpp"
#include "Shader.hpp"
#include <vector>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

struct Renderer {
    Renderer();
    bool draw(const Scene& scene);

private:
    GLFWwindow* window;
    int screen_width = 800;
    int screen_height = 600;
    Camera camera;
    float lastFrameTime = 0.0f;
    float deltaTime = 0.0f;
    Shader bboxShader, particleShader;
    GLuint bbox_VAO, bbox_VBO, bbox_EBO;
    GLuint particle_VAO, particle_VBO, particle_instanceVBO;
    vector<vec3> sphere_vertices;
    bool mouseContained = false;

    void framebufferSizeCallback(GLFWwindow* window, int width, int height);
    void processInput(GLFWwindow* window);
    void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
    void mouseCallback(GLFWwindow* window, double xpos, double ypos);
    void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
    void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);

    void drawBbox(const Scene& scene);
    void drawParticles(const Scene& scene);
};
