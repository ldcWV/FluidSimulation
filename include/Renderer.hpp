#pragma once
#include "Scene.hpp"
#include <GLFW/glfw3.h>

struct Renderer {
    Renderer(GLFWwindow* window);
    void draw(const Scene& scene);

private:
    GLFWwindow* window;
    int width, height;
};
