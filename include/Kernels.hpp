#pragma once
#include "glm/glm.hpp"

namespace Kernels {
    double poly6(const glm::dvec3& r, const double h);
    glm::dvec3 gradSpiky(const glm::dvec3& r, const double h);
};
