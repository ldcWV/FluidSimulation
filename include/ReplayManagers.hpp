#pragma once
#include <string>
#include <Scene.hpp>

struct ReplayWriter {
    ReplayWriter(std::string fname);
    ~ReplayWriter();
    void write_scene(const Scene& scene);

    FILE* file;
};

struct ReplayReader {
    ReplayReader(std::string fname);
    ~ReplayReader();
    int read_scene(Scene* scene);
    void reset();

    FILE* file;
};
