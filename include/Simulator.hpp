#pragma once
#include "Scene.hpp"

struct Simulator {
    virtual void update(double elapsed, Scene& scene) = 0;
};

struct SequentialSimulator : Simulator {
    SequentialSimulator(const Scene& scene);
    ~SequentialSimulator();

    void update(double elapsed, Scene& scene) override;

private:
    size_t grid_width, grid_height, grid_length;
    int* grid;
    int* grid_cell_counts;
    int* neighbors;
    int* neighbor_counts;
    glm::dvec3 bbox_mins;
    glm::dvec3 bbox_maxs;

    glm::ivec3 get_cell_coords(glm::dvec3 pos);
    int get_cell_idx(glm::ivec3 coords);
    void recompute_grid(const Scene& scene);
    void recompute_neighbors(const Scene& scene);
};

struct ParallelSimulator : Simulator {
    ParallelSimulator(const Scene& scene);
    ~ParallelSimulator();

    void update(double elapsed, Scene& scene) override;
};
