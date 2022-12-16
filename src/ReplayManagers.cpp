#include "ReplayManagers.hpp"
#include <iostream>

using namespace std;

ReplayWriter::ReplayWriter(string fname) {
    file = fopen(fname.c_str(), "wb");
    if (file == nullptr) {
        cout << "Failed to read binary scene file!" << endl;
        return;
    }
}

ReplayWriter::~ReplayWriter() {
    fclose(file);
}

void ReplayWriter::write_scene(const Scene& scene) {
    int sz = scene.particles.size();
    fwrite(&sz, sizeof(int), 1, file);

    fwrite(&scene.bbox_mins[0], sizeof(glm::dvec3), 1, file);
    fwrite(&scene.bbox_maxs[0], sizeof(glm::dvec3), 1, file);

    fwrite(&scene.particles[0], sizeof(Particle), scene.particles.size(), file);
}

ReplayReader::ReplayReader(string fname) {
    file = fopen(fname.c_str(), "rb");
    if (file == nullptr) {
        cout << "Failed to read binary scene file!" << endl;
        return;
    }
}

ReplayReader::~ReplayReader() {
    fclose(file);
}

int ReplayReader::read_scene(Scene* scene) {
    int sz;
    if (fread(&sz, sizeof(int), 1, file) < 1) return 1;

    if (fread(&(scene->bbox_mins[0]), sizeof(glm::dvec3), 1, file) < 1) return 2;
    if (fread(&(scene->bbox_maxs[0]), sizeof(glm::dvec3), 1, file) < 1) return 3;

    scene->particles.resize(sz);
    if (fread(&(scene->particles[0]), sizeof(Particle), sz, file) < sz) return 4;

    return 0;
}

void ReplayReader::reset() {
    rewind(file);
}
