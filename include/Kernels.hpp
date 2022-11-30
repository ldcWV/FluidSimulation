#pragma once
#include "glm/glm.hpp"
#include <math.h>

namespace Kernels {
    double poly6(const glm::dvec3& r, const double h);

    glm::dvec3 gradSpiky(const glm::dvec3& r, const double h) {
        double r_mag = glm::length(r);
        if (r_mag > h) return glm::dvec3{0, 0, 0};
        return -45 / (Constants::pi * pow(h, 6)) * pow(h - r_mag, 2) * glm::normalize(r);
    }
};
