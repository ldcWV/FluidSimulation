#pragma once
#include "Scene.hpp"
#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include "Camera.hpp"
#include "Shader.hpp"

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
    Shader particleShader;
    Shader bboxShader;
    GLuint bbox_VAO;

    void framebufferSizeCallback(GLFWwindow* window, int width, int height);
    void processInput(GLFWwindow* window);
    void mouseCallback(GLFWwindow* window, double xpos, double ypos);
    void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);

    void drawBbox(const Scene& scene);
};
