#include "Simulator.hpp"
#include "Constants.hpp"
#include "Kernels.hpp"
#include <iostream>

using namespace glm;

SequentialSimulator::SequentialSimulator(const Scene& scene) {
    this->bbox_mins = scene.bbox_mins;
    this->bbox_maxs = scene.bbox_maxs;

    this->grid_width = (bbox_maxs.x - bbox_mins.x) / Constants::h + 1;
    this->grid_height = (bbox_maxs.y - bbox_mins.y) / Constants::h + 1;
    this->grid_length = (bbox_maxs.z - bbox_mins.z) / Constants::h + 1;

    this->grid = (int*)calloc(grid_width * grid_height * grid_length * Constants::max_particles_per_cell, sizeof(int));
    this->grid_cell_counts = (int*)calloc(grid_width * grid_height * grid_length, sizeof(int));
    this->neighbors = (int*)calloc(scene.particles.size()*Constants::max_neighbors, sizeof(int));
    this->neighbor_counts = (int*)calloc(scene.particles.size(), sizeof(int));
}

SequentialSimulator::~SequentialSimulator() {
    free(grid);
    free(grid_cell_counts);
    free(neighbors);
    free(neighbor_counts);
}

ivec3 SequentialSimulator::get_cell_coords(dvec3 pos) {
    return ivec3{
        (pos.x - bbox_mins.x) / Constants::h,
        (pos.y - bbox_mins.y) / Constants::h,
        (pos.z - bbox_mins.z) / Constants::h
    };
}

int SequentialSimulator::get_cell_idx(ivec3 coords) {
    return coords.x * grid_height*grid_length + coords.y * grid_length + coords.z;
}

void SequentialSimulator::recompute_grid(const Scene& scene) {
    for (int i = 0; i < grid_width*grid_height*grid_length; i++) {
        grid_cell_counts[i] = 0;
    }
    for (auto& p : scene.particles) {
        ivec3 coords = get_cell_coords(p.new_pos);
        int idx = get_cell_idx(coords);
        if (grid_cell_counts[idx] == Constants::max_particles_per_cell) continue;
        grid[idx * Constants::max_particles_per_cell + grid_cell_counts[idx]++] = p.id;
    }
}

void SequentialSimulator::recompute_neighbors(const Scene& scene) {
    for (auto& p : scene.particles) {
        ivec3 our_coords = get_cell_coords(p.new_pos);

        // Need to check the 27 surrounding cells
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    ivec3 coords = our_coords + ivec3{dx, dy, dz};
                    if (coords.x < 0 || coords.x >= grid_width ||
                        coords.y < 0 || coords.y >= grid_height ||
                        coords.z < 0 || coords.z >= grid_length) continue;
                    
                    // Iterate through all particles here (including ourselves!) and add it to neighbors
                    int idx = get_cell_idx(coords);
                    for (int i = 0; i < grid_cell_counts[idx]; i++) {
                        int other_id = grid[idx*Constants::max_particles_per_cell + i];
                        if (neighbor_counts[p.id] < Constants::max_neighbors) {
                            neighbors[p.id*Constants::max_neighbors + neighbor_counts[p.id]++] = other_id;
                        }
                    }
                }
            }
        }
    }
}

void SequentialSimulator::update(double elapsed, Scene& scene) {
    for (auto& p : scene.particles) {
        // Apply forces
        p.vel += elapsed * Constants::g;

        // Predict new position
        p.new_pos = p.pos + elapsed * p.vel;
        // Collisions with box
        for (int i = 0; i < 3; i++) {
            if (p.new_pos[i] < bbox_mins[i] + Constants::radius) {
                p.vel[i] = 0;
                p.new_pos[i] = bbox_mins[i] + Constants::radius;
            } else if (p.new_pos[i] > bbox_maxs[i] - Constants::radius) {
                p.vel[i] = 0;
                p.new_pos[i] = bbox_maxs[i] - Constants::radius;
            }
        }
    }

    recompute_grid(scene);
    recompute_neighbors(scene);

    for (int i = 0; i < Constants::solver_iterations; i++) {
        // Calculate lambda_i

        // Calculate delta_p and perform collision detection and response

        // Update new_pos
    }

    for (auto& p : scene.particles) {
        p.vel = 1.0 / elapsed * (p.new_pos - p.pos);
        // TODO: apply vorticity confinement and XSPH viscosity
        p.pos = p.new_pos;
    }
}