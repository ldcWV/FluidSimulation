#pragma once

#include <glad/gl.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

enum CameraMovement {
    FORWARD,
    BACKWARD,
    LEFT,
    RIGHT,
    UP,
    DOWN
};

class Camera {
public:
    glm::vec3 Position;
    glm::vec3 Front;
    glm::vec3 Right;
    glm::vec3 WorldUp;
    glm::vec3 Up;
    float Yaw;
    float Pitch;
    float MovementSpeed;
    float MouseSensitivity;
    float ScrollSensitivity;
    float FOV;

    Camera(glm::vec3 position = glm::vec3(0.f, 0.f, 0.f), glm::vec3 up = glm::vec3(0.f, 1.f, 0.f), float yaw = YAW, float pitch = PITCH);
    Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw = YAW, float pitch = PITCH);
    glm::mat4 getViewMatrix();
    void processKeyboard(CameraMovement direction, float elapsed);
    void processMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch = true);
    void processMouseScroll(float yoffset);

private:
    // defaults
    static constexpr float YAW = -90.f;
    static constexpr float PITCH = 0.0f;
    static constexpr float SPEED = 5.f;
    static constexpr float SENSITIVITY = 0.05f;
    static constexpr float SCROLL_SENSITIVITY = 4.f;
    static constexpr float defaultFOV = 60.0f;

    void updateCameraVectors();
};
