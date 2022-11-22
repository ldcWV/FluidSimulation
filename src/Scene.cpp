#include "Scene.hpp"
#include <fstream>

using namespace std;

Scene::Scene() {
}

Scene::Scene(string filename) {
    ifstream fin(filename.c_str());

    int num_particles;
    fin >> num_particles;

    for (int i = 0; i < 3; i++) fin >> bbox_mins[i];
    for (int i = 0; i < 3; i++) fin >> bbox_maxs[i];
    
    for (int i = 0; i < num_particles; i++) {
        Particle p;
        fin >> p.id;
        for (int j = 0; j < 3; j++) fin >> p.old_pos[j];
        for (int j = 0; j < 3; j++) fin >> p.pos[j];
        for (int j = 0; j < 3; j++) fin >> p.vel[j];
        particles.push_back(p);
    }
}

void Scene::save(string filename) {
    ofstream fout(filename.c_str());

    fout << particles.size() << "\n";

    for (int i = 0; i < 3; i++) fout << bbox_mins[i];
    for (int i = 0; i < 3; i++) fout << bbox_maxs[i];

    for (Particle p : particles) {
        fout << p.id;
        for (int j = 0; j < 3; j++) fout << p.old_pos[j];
        for (int j = 0; j < 3; j++) fout << p.pos[j];
        for (int j = 0; j < 3; j++) fout << p.vel[j];
    }
}
