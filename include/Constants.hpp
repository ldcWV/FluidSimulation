#pragma once
#include "glm/glm.hpp"

namespace Constants {
    const glm::dvec3 g = {0, -9.8, 0};
    const double dt = 1.0 / 60.0 * 0.2;
    const double eps = 1e-5;
    const double mass = 1.0 / 30.0;
    const double radius = 0.02;
    const double h = 0.1; // not the same thing as radius 
    const size_t max_particles_per_cell = 100;
    const size_t max_neighbors = 100;
    const int solver_iterations = 3;
    const double pi = 3.14159265358979323846;
    const double rest_density = 60;
    const double corr_q = 0.3; 
    const double corr_k = 0.0001;
    const int corr_n = 4;
    const float xsph_c = 0.001f;
};
