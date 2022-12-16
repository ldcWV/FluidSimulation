#include "Scene.hpp"
#include <fstream>
#include <iostream>

using namespace std;

Scene::Scene() {
}

Scene::Scene(const Scene& other) {
    this->bbox_mins = other.bbox_mins;
    this->bbox_maxs = other.bbox_maxs;
    this->particles = other.particles; 
}

Scene::Scene(string filename) {
    ifstream fin(filename.c_str());
    if (fin.fail()) {
        cout << "Scene load failed!" << endl;
        return;
    }

    int num_particles;
    fin >> num_particles;

    for (int i = 0; i < 3; i++) fin >> bbox_mins[i];
    for (int i = 0; i < 3; i++) fin >> bbox_maxs[i];
    
    particles.clear();
    for (int i = 0; i < num_particles; i++) {
        Particle p;
        fin >> p.id;
        for (int j = 0; j < 3; j++) fin >> p.pos[j];
        for (int j = 0; j < 3; j++) fin >> p.new_pos[j];
        for (int j = 0; j < 3; j++) fin >> p.vel[j];
        particles.push_back(p);
    }
}
