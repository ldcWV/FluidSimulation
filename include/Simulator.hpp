#pragma once
#include "Scene.hpp"
#include "Particle.hpp"

struct Simulator {
    virtual void update(double elapsed, Scene& scene) = 0;
};

struct SequentialSimulator : Simulator {
    SequentialSimulator(const Scene& scene);
    ~SequentialSimulator();

    void update(double elapsed, Scene& scene) override;

private:
    size_t grid_width, grid_height, grid_length;
    int* grid_sizes;
    int* grid_offsets;
    int* grid;
    int* neighbor_sizes;
    int* neighbor_offsets;
    static constexpr int MAX_NEIGHBORS = 5000000; // sum of # neighbors over every particle
    int neighbors[MAX_NEIGHBORS];
    double* lambdas;
    double* densities;
    glm::dvec3* delta_pos;
    glm::dvec3* delta_vel;
    glm::dvec3 bbox_mins;
    glm::dvec3 bbox_maxs;

    glm::ivec3 get_cell_coords(glm::dvec3 pos);
    int get_cell_idx(glm::ivec3 coords);
    void recompute_grid(Scene& scene);
    void recompute_neighbors(Scene& scene);
    double compute_density(Scene& scene, int particle_id);
    double compute_constraint(Scene& scene, int particle_id);
    glm::dvec3 compute_grad_constraint(Scene& scene, int constraint_id, int grad_id);
};

struct ParallelSimulator : Simulator {
    ParallelSimulator(const Scene& scene);
    ~ParallelSimulator();
    void update(double elapsed, Scene& scene) override;

private: 
    // MAKE SURE THIS STUFF ISN'T OVERLOADED??? 

    // HOST memory 
    int _n;
    size_t _blocks, _threads;
    size_t _total_cells; 

    // DEVICE memory 
    size_t n, total_cells; 
    glm::dvec3 bbox_mins, bbox_maxs;
    size_t grid_width, grid_height, grid_length;
    double* lambdas;
    double* densities; 
    Particle* particles;
    int* grid;
    int* neighbors;
    int *bins, *prefix_bins, *neighbor_starts, *neighbor_sizes, *grid_starts;
    glm::dvec3* delta_pos, *delta_vel;

    // HOST functions 
    void reset(); 
    void recompute_neighbors();
    void recompute_grid();
    void simulate(); 

    // HOST kernel wrappers
    void compute_lambdas();
    void compute_delta_positions();
    void compute_densities(); 
    void compute_velocities(double elapsed);
    void compute_velocities_and_positions(double elapsed);
    void xsph_viscosity();
    void update_velocities();
    void update_positions();
    void update_collisions();
};