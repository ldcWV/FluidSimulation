#pragma once
#include "glm/glm.hpp"

namespace Constants {
    const glm::dvec3 g = {0, -9.8, 0};
    const double dt = 1.0 / 6000.0;
    const double eps = 1e5;
    const double mass = 1.0 / 30.0;
    const double radius = 0.1;
    const double h = 1.0;
    const size_t max_particles_per_cell = 1000;
    const size_t max_neighbors = 100;
    const int solver_iterations = 3;
    const double pi = 3.14159265358979323846;
    const double rest_density = 1000;
};
