#include "Renderer.hpp"
#include <GLFW/glfw3.h>
#include <math.h>
#include <iostream>
#include "Constants.hpp"

void DrawCircle(float cx, float cy, float r, int num_segments) {
    glBegin(GL_LINE_LOOP);
    for(int ii = 0; ii < num_segments; ii++)
    {
        float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle

        float x = r * cosf(theta);//calculate the x component
        float y = r * sinf(theta);//calculate the y component

        glVertex2f(x + cx, y + cy);//output vertex

    }
    glEnd();
}

void DrawLine(float x0, float y0, float x1, float y1) {
    glBegin(GL_LINES);
    glVertex2f(x0, y0);
    glVertex2f(x1, y1);
    glEnd();
}

Renderer::Renderer(GLFWwindow* window) : window(window) {
    glfwGetWindowSize(window, &width, &height);
    glfwMakeContextCurrent(window);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
}

void Renderer::draw(const Scene& scene) {
    // Set up camera position
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    float aspect = 1.f * width / height;
    float min_width = 1.2f*(scene.bbox_maxs.x - scene.bbox_mins.x);
    float min_height = 1.2f*(scene.bbox_maxs.y - scene.bbox_mins.y);
    float draw_width = max(min_width, min_height * aspect);
    float draw_height = draw_width / aspect;
    glm::dvec3 mid = (scene.bbox_mins + scene.bbox_maxs) / 2.;
    glOrtho(
        mid.x - draw_width/2, mid.x + draw_width/2,
        mid.y - draw_height/2, mid.y + draw_height/2,
        -1.0, 1.0
    );

    // Clear canvas
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw bounding box
    glColor3f(0.0, 0.0, 0.0);
    DrawLine(scene.bbox_mins.x, scene.bbox_mins.y, scene.bbox_maxs.x, scene.bbox_mins.y);
    DrawLine(scene.bbox_maxs.x, scene.bbox_mins.y, scene.bbox_maxs.x, scene.bbox_maxs.y);
    DrawLine(scene.bbox_maxs.x, scene.bbox_maxs.y, scene.bbox_mins.x, scene.bbox_maxs.y);
    DrawLine(scene.bbox_mins.x, scene.bbox_maxs.y, scene.bbox_mins.x, scene.bbox_mins.y);

    // Draw particles
    glColor3f(1.0, 0.0, 0.0);
    for (auto particle : scene.particles) {
        DrawCircle(particle.pos.x, particle.pos.y, Constants::radius, 64);
    }
}
