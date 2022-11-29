#pragma once

#include "glm/glm.hpp"

struct Particle {
    int id;
    glm::dvec3 old_pos;
    glm::dvec3 pos;
    glm::dvec3 vel;
};
