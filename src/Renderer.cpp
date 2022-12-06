#include "Renderer.hpp"
#include <math.h>
#include <iostream>
#include "Constants.hpp"
#include "Shader.hpp"
#include "Camera.hpp"

using namespace glm;
using namespace std;

void Renderer::framebufferSizeCallback(GLFWwindow* window, int width, int height) {
    screen_width = width;
    screen_height = height;
    glViewport(0, 0, width, height);
}

void Renderer::processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, true);
    }

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camera.processKeyboard(CameraMovement::FORWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camera.processKeyboard(CameraMovement::BACKWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camera.processKeyboard(CameraMovement::LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camera.processKeyboard(CameraMovement::RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
        camera.processKeyboard(CameraMovement::UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS) {
        camera.processKeyboard(CameraMovement::DOWN, deltaTime);
    }
}

void Renderer::mouseCallback(GLFWwindow* window, double xpos, double ypos) {
    static bool firstMouse = true;
    static double lastX = 0, lastY = 0;

    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = float(xpos - lastX);
    float yoffset = float(lastY - ypos);
    camera.processMouseMovement(xoffset, yoffset, true);

    lastX = xpos;
    lastY = ypos;
}

void Renderer::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    camera.processMouseScroll((float)yoffset);
}

Renderer::Renderer() {
    // Initialize window
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(screen_width, screen_height, "Epic Fluid Simulation", NULL, NULL);
    if (window == nullptr) {
        cout << "Failed to create GLFW window" << endl;
        return;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoaderLoadGL()) {
        cout << "Failed to initialize GLAD" << endl;
        return;
    }

    glViewport(0, 0, screen_width, screen_height);

    // Set up various callbacks
    glfwSetWindowUserPointer(window, this);
    auto framebufferSizeCallbackLambda = [](GLFWwindow* window, int width, int height) {
        static_cast<Renderer*>(glfwGetWindowUserPointer(window))->framebufferSizeCallback(window, width, height);
    };
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallbackLambda); 

    auto mouseCallbackLambda = [](GLFWwindow* window, double xpos, double ypos) {
        static_cast<Renderer*>(glfwGetWindowUserPointer(window))->mouseCallback(window, xpos, ypos);
    };
    glfwSetCursorPosCallback(window, mouseCallbackLambda);

    auto scrollCallbackLambda = [](GLFWwindow* window, double xoffset, double yoffset) {
        static_cast<Renderer*>(glfwGetWindowUserPointer(window))->scrollCallback(window, xoffset, yoffset);
    };
    glfwSetScrollCallback(window, scrollCallbackLambda);

    // Disappear cursor and keep it in the window
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    camera = Camera(vec3(0.f, 0.f, 20.f));

    glEnable(GL_DEPTH_TEST);

    std::string vPath = std::string(SRC_DIR) + "ParticleShader.vert";
    std::string fPath = std::string(SRC_DIR) + "ParticleShader.frag";
    particleShader = Shader(vPath.c_str(), fPath.c_str());

    vPath = std::string(SRC_DIR) + "BboxShader.vert";
    fPath = std::string(SRC_DIR) + "BboxShader.frag";
    bboxShader = Shader(vPath.c_str(), fPath.c_str());

    // Prepare bbox_VAO
    glGenVertexArrays(1, &bbox_VAO);
    GLuint bbox_VBO, bbox_EBO;
    glGenBuffers(1, &bbox_VBO);
    glGenBuffers(1, &bbox_EBO);

    glBindVertexArray(bbox_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, bbox_VBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bbox_EBO);

    unsigned int indices[] = {
        0, 1, 1, 2, 2, 3, 3, 0, // bottom face
        0, 4, 1, 5, 2, 6, 3, 7, // sides
        4, 5, 5, 6, 6, 7, 7, 4 // top face
    };

    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);
}

bool Renderer::draw(const Scene& scene) {
    if (glfwWindowShouldClose(window)) {
        cout << "Closing glfw window" << endl;
        return false;
    }

    processInput(window);

    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // glm::mat4 model(1.f);
    // glm::mat4 view = camera.getViewMatrix();
    // glm::mat4 projection = glm::perspective(glm::radians(camera.FOV), float(screen_width/screen_height), 0.1f, 100.f);

    drawBbox(scene);

    glfwSwapBuffers(window);
    glfwPollEvents();

    float currentFrameTime = glfwGetTime();
    deltaTime = currentFrameTime - lastFrameTime;
    lastFrameTime = currentFrameTime;

    return true;
}

void Renderer::drawBbox(const Scene& scene) {
    vec3 mins = scene.bbox_mins;
    vec3 maxs = scene.bbox_maxs;

    float vertices[] = {
        mins.x, mins.y, mins.z,
        maxs.x, mins.y, mins.z,
        maxs.x, mins.y, maxs.z,
        mins.x, mins.y, maxs.z,
        mins.x, maxs.y, mins.z,
        maxs.x, maxs.y, mins.z,
        maxs.x, maxs.y, maxs.z,
        mins.x, maxs.y, maxs.z
    };

    glBindVertexArray(bbox_VAO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);

    glm::mat4 view = camera.getViewMatrix();
    glm::mat4 projection = glm::perspective(glm::radians(camera.FOV), float(screen_width/screen_height), 0.1f, 100.f);
    glm::vec4 test = projection * view * glm::vec4(vertices[0], vertices[1], vertices[2], 1.f);
    bboxShader.use();
    bboxShader.setMat4("view", view);
    bboxShader.setMat4("projection", projection);

    glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
}
