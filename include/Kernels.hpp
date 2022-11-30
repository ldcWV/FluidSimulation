#pragma once
#include "glm/glm.hpp"
#include <math.h>

namespace Kernels {
    double poly6(const glm::dvec3& r, const double h) {
        double r_mag = glm::length(r);
        if (r_mag > h) 0;
        return 315.0 / (64 * Constants::pi * pow(h, 9)) * pow(h*h - r_mag*r_mag, 3);
    }

    glm::dvec3 gradSpiky(const glm::dvec3& r, const double h) {
        double r_mag = glm::length(r);
        if (r_mag > h) return glm::dvec3{0, 0, 0};
        return -45 / (Constants::pi * pow(h, 6)) * pow(h - r_mag, 2) * glm::normalize(r);
    }
};
