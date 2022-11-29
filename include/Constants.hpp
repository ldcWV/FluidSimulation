#pragma once
#include "glm/glm.hpp"

namespace Constants {
    const glm::dvec3 g = {0, -9.8, 0};
    const double dt = 1.0 / 60.0;
    const double eps = 1e5;
    const double mass = 1.0 / 30.0;
    const double radius = 0.1;
};