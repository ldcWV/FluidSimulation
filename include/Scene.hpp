#pragma once
#include <vector>
#include "Particle.hpp"
#include <string>

using namespace std;

struct Scene {
    glm::dvec3 bbox_mins = {-5, -5, -5};
    glm::dvec3 bbox_maxs = { 5,  5,  5};
    vector<Particle> particles;

    Scene();
    Scene(string fname, bool binary);
    int load(string fname, bool binary);
    void save(string fname);
};
