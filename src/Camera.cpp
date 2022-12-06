#include "Camera.hpp"

// Source: learnopengl.com

Camera::Camera(glm::vec3 position, glm::vec3 up, float yaw, float pitch) {
    Position = position;
    WorldUp = up;
    Yaw = yaw;
    Pitch = pitch;
    MovementSpeed = SPEED;
    MouseSensitivity = SENSITIVITY;
    ScrollSensitivity = SCROLL_SENSITIVITY;
    FOV = defaultFOV;
    updateCameraVectors();
}

Camera::Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) {
    Position = glm::vec3(posX, posY, posZ);
    WorldUp = glm::vec3(upX, upY, upZ);
    Yaw = yaw;
    Pitch = pitch;
    MovementSpeed = SPEED;
    MouseSensitivity = SENSITIVITY;
    ScrollSensitivity = SCROLL_SENSITIVITY;
    FOV = defaultFOV;
    updateCameraVectors();
}

glm::mat4 Camera::getViewMatrix() {
    return glm::lookAt(Position, Position + Front, WorldUp);
}

void Camera::processKeyboard(CameraMovement direction, float elapsed) {
    float amt = MovementSpeed * elapsed;
    if (direction == CameraMovement::FORWARD) {
        Position += amt * Front;
    } else if (direction == CameraMovement::BACKWARD) {
        Position -= amt * Front;
    } else if (direction == CameraMovement::LEFT) {
        Position -= amt * Right;
    } else if (direction == CameraMovement::RIGHT) {
        Position += amt * Right;
    } else if (direction == CameraMovement::UP) {
        Position += amt * WorldUp;
    } else if (direction == CameraMovement::DOWN) {
        Position -= amt * WorldUp;
    }
}

void Camera::processMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch) {
    Yaw += xoffset * MouseSensitivity;
    Pitch += yoffset * MouseSensitivity;

    if (constrainPitch) {
        Pitch = fmax(Pitch, -89.f);
        Pitch = fmin(Pitch, 89.f);
    }

    updateCameraVectors();
}

void Camera::processMouseScroll(float yoffset) {
    FOV -= yoffset * ScrollSensitivity;
    FOV = fmax(FOV, 1.f);
    FOV = fmin(FOV, defaultFOV);
}

void Camera::updateCameraVectors() {
    Front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
    Front.y = sin(glm::radians(Pitch));
    Front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
    Front = glm::normalize(Front);
    Right = glm::normalize(glm::cross(Front, WorldUp));
    Up = glm::cross(Right, Front);
}
