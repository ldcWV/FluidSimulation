#include "Scene.hpp"
#include <fstream>
#include <iostream>

using namespace std;

Scene::Scene() {
}

Scene::Scene(string filename, bool binary) {
    load(filename, binary);
}

int Scene::load(string filename, bool binary) {
    if (binary) {
        FILE* file = fopen(filename.c_str(), "rb");
        if (file == nullptr) {
            return -1;
        }

        int num_particles;
        fread(&num_particles, sizeof(int), 1, file);

        fread(&bbox_mins[0], sizeof(glm::dvec3), 1, file);
        fread(&bbox_maxs[0], sizeof(glm::dvec3), 1, file);

        particles.resize(num_particles);
        fread(&particles[0], sizeof(Particle), particles.size(), file);
        fclose(file);
    } else {
        ifstream fin(filename.c_str());
        if (fin.fail()) {
            return -1;
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
    return 0;
}

void Scene::save(string filename) {
    FILE* file = fopen(filename.c_str(), "wb");
    if (file == nullptr) {
        cout << "Failed to read binary scene file!" << endl;
        return;
    }
    int sz = particles.size();
    fwrite(&sz, sizeof(int), 1, file);

    fwrite(&bbox_mins[0], sizeof(glm::dvec3), 1, file);
    fwrite(&bbox_maxs[0], sizeof(glm::dvec3), 1, file);

    fwrite(&particles[0], sizeof(Particle), particles.size(), file);
    fclose(file);
}
