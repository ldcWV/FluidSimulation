#pragma once

#include "glm/glm.hpp"

struct Particle {
    int id;
    glm::dvec3 pos;
    glm::dvec3 new_pos;
    glm::dvec3 vel;
    int mark; // stores info for debugging
};
