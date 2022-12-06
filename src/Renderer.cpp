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

    camera = Camera(vec3(0.f, 0.f, 5.f));

    glEnable(GL_DEPTH_TEST);

    float vertices[] = {
        -0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.0f,
         0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
         0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
        -0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 1.0f,
        -0.5f, -0.5f,  0.5f, 1.0f, 0.0f, 0.0f,
         0.5f, -0.5f,  0.5f, 1.0f, 0.0f, 1.0f,
         0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
        -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f
    };

    unsigned int indices[] = {
        0, 1, 2, // bottom
        0, 2, 3,
        0, 1, 5, // front
        0, 4, 5,
        1, 5, 6, // right
        1, 2, 6,
        2, 3, 6, // back
        3, 6, 7,
        0, 3, 7, // left
        0, 4, 7,
        4, 5, 6, // top
        4, 6, 7
    };

    GLuint VAO;
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    GLuint VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    GLuint EBO;
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(3 * sizeof(float)));

    std::string vShaderPath = std::string(SRC_DIR) + "shader.vert";
    std::string fShaderPath = std::string(SRC_DIR) + "shader.frag";
    shader = Shader(vShaderPath.c_str(), fShaderPath.c_str());
    shader.use();
}

bool Renderer::draw(const Scene& scene) {
    if (glfwWindowShouldClose(window)) {
        cout << "Closing glfw window" << endl;
        return false;
    }

    processInput(window);

    glm::vec3 cubePositions[] = {
        glm::vec3( 0.0f,  0.0f,  0.0f), 
        glm::vec3( 2.0f,  5.0f, -15.0f), 
        glm::vec3(-1.5f, -2.2f, -2.5f),  
        glm::vec3(-3.8f, -2.0f, -12.3f),  
        glm::vec3( 2.4f, -0.4f, -3.5f),  
        glm::vec3(-1.7f,  3.0f, -7.5f),  
        glm::vec3( 1.3f, -2.0f, -2.5f),  
        glm::vec3( 1.5f,  2.0f, -2.5f), 
        glm::vec3( 1.5f,  0.2f, -1.5f), 
        glm::vec3(-1.3f,  1.0f, -1.5f)  
    };
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glm::mat4 view = camera.getViewMatrix();
    glm::mat4 projection = glm::perspective(glm::radians(camera.FOV), float(screen_width/screen_height), 0.1f, 100.f);

    shader.setMat4("view", view);
    shader.setMat4("projection", projection);

    for (int i = 0; i < 10; i++) {
        glm::mat4 model(1.f);
        model = glm::translate(model, cubePositions[i]);
        float angle = glm::radians(20.f * i) + glm::radians(50.f) * glfwGetTime();
        model = glm::rotate(model, angle, glm::vec3(1.0f, 0.3f, 0.5f));
        shader.setMat4("model", model);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    }

    glfwSwapBuffers(window);
    glfwPollEvents();

    float currentFrameTime = glfwGetTime();
    deltaTime = currentFrameTime - lastFrameTime;
    lastFrameTime = currentFrameTime;

    return true;
}
