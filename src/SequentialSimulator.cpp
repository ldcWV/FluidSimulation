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

    this->grid_sizes = (int*)calloc(grid_width * grid_height * grid_length, sizeof(int));
    this->grid_offsets = (int*)calloc(grid_width * grid_height * grid_length, sizeof(int));
    this->grid = (int*)calloc(scene.particles.size(), sizeof(int));

    this->neighbor_sizes = (int*)calloc(scene.particles.size(), sizeof(int));
    this->neighbor_offsets = (int*)calloc(scene.particles.size(), sizeof(int));

    this->lambdas = (double*)calloc(scene.particles.size(), sizeof(double));
    this->densities = (double*)calloc(scene.particles.size(), sizeof(double));
    this->delta_pos = (dvec3*)calloc(scene.particles.size(), sizeof(dvec3));
    this->delta_vel = (dvec3*)calloc(scene.particles.size(), sizeof(dvec3));
}

SequentialSimulator::~SequentialSimulator() {
    free(lambdas);
    free(densities);
    free(delta_pos);
    free(delta_vel);
    free(grid_sizes);
    free(grid_offsets);
    free(grid);
    free(neighbor_sizes);
    free(neighbor_offsets);
}

ivec3 SequentialSimulator::get_cell_coords(dvec3 pos) {
    ivec3 res{
        (pos.x - bbox_mins.x) / Constants::h,
        (pos.y - bbox_mins.y) / Constants::h,
        (pos.z - bbox_mins.z) / Constants::h
    };
    return clamp(res, ivec3{0, 0, 0}, ivec3{grid_width-1, grid_height-1, grid_length-1});
}

int SequentialSimulator::get_cell_idx(ivec3 coords) {
    return coords.x * grid_height*grid_length + coords.y * grid_length + coords.z;
}

void SequentialSimulator::recompute_grid(Scene& scene) {
    // Compute sizes and offsets of each grid cell
    for (size_t i = 0; i < grid_width*grid_height*grid_length; i++) {
        grid_sizes[i] = 0;
    }
    for (auto& p : scene.particles) {
        ivec3 coords = get_cell_coords(p.new_pos);
        int idx = get_cell_idx(coords);
        grid_sizes[idx]++;
    }
    grid_offsets[0] = 0;
    for (size_t i = 1; i < grid_width*grid_height*grid_length; i++) {
        grid_offsets[i] = grid_sizes[i-1] + grid_offsets[i-1];
    }

    // Fill in the grid cells using the computed offsets
    for (size_t i = 0; i < grid_width*grid_height*grid_length; i++) {
        grid_sizes[i] = 0;
    }
    for (auto& p : scene.particles) {
        ivec3 coords = get_cell_coords(p.new_pos);
        int idx = get_cell_idx(coords);
        grid[grid_offsets[idx] + grid_sizes[idx]] = p.id;
        grid_sizes[idx]++;
    }
}

void SequentialSimulator::recompute_neighbors(Scene& scene) {
    // Compute sizes and contents of each neighbor list
    for (size_t i = 0; i < scene.particles.size(); i++) {
        neighbor_sizes[i] = 0;
    }
    int n_idx = 0;
    for (size_t i = 0; i < scene.particles.size(); i++) {
        Particle& p = scene.particles[i];
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
                    for (int j = grid_offsets[idx]; j < grid_offsets[idx] + grid_sizes[idx]; j++) {
                        int other_id = grid[j];
                        if (n_idx == MAX_NEIGHBORS) {
                            cout << "Warning: ran out of space for neighbors (Constants::MAX_NEIGHBORS should be increased)." << endl;
                            continue;
                        }
                        neighbors[n_idx] = other_id;
                        n_idx++;
                        neighbor_sizes[i]++;
                    }
                }
            }
        }
    }

    // Compute offsets of each neighbor list
    neighbor_offsets[0] = 0;
    for (size_t i = 1; i < scene.particles.size(); i++) {
        neighbor_offsets[i] = neighbor_sizes[i-1] + neighbor_offsets[i-1];
    }
}

double SequentialSimulator::compute_density(Scene& scene, int particle_id) {
    const Particle& target = scene.particles[particle_id];
    double density = 0.0;
    for (int i = neighbor_offsets[particle_id]; i < neighbor_offsets[particle_id] + neighbor_sizes[particle_id]; i++) {
        int neighbor_id = neighbors[i];
        const Particle& p = scene.particles[neighbor_id];
        density += Constants::mass * Kernels::poly6(target.new_pos - p.new_pos, Constants::h);
    }
    return density;
}

double SequentialSimulator::compute_constraint(Scene& scene, int particle_id) {
    return compute_density(scene, particle_id) / Constants::rest_density - 1.0;
}

// Assumes that grad_id is the id of a neighbor of particle_id
dvec3 SequentialSimulator::compute_grad_constraint(Scene& scene, int constraint_id, int grad_id) {
    const Particle& constraint_particle = scene.particles[constraint_id];
    const dvec3& constraint_pos = constraint_particle.new_pos;
    if (constraint_id == grad_id) {
        dvec3 res{0.0, 0.0, 0.0};
        for (int i = neighbor_offsets[constraint_id]; i < neighbor_offsets[constraint_id] + neighbor_sizes[constraint_id]; i++) {
            int neighbor_id = neighbors[i];
            const dvec3& neighbor_pos = scene.particles[neighbor_id].new_pos;
            res += Constants::mass * Kernels::gradSpiky(constraint_pos - neighbor_pos, Constants::h);
        }
        return res / Constants::rest_density;
    } else {
        return -Constants::mass * Kernels::gradSpiky(constraint_pos - scene.particles[grad_id].new_pos, Constants::h) / Constants::rest_density;
    }
}

void SequentialSimulator::update(double elapsed, Scene& scene) {
    for (auto& p : scene.particles) {
        // Apply forces
        p.vel += elapsed * Constants::g;
        // Predict new position
        p.new_pos = p.pos + elapsed * p.vel;
    }

    recompute_grid(scene);
    recompute_neighbors(scene);
    
    for (int iter = 0; iter < Constants::solver_iterations; iter++) {
        // Calculate lambda_i
        for (int i = 0; i < scene.particles.size(); ++i) {
            double numerator = compute_constraint(scene, i);

            double denominator = 0.0;
            for (int j = neighbor_offsets[i]; j < neighbor_offsets[i] + neighbor_sizes[i]; j++) {
                int neighbor_id = neighbors[j];
                auto grad = compute_grad_constraint(scene, i, neighbor_id);
                denominator += dot(grad, grad) / Constants::mass;
            }
            denominator += Constants::eps;

            lambdas[i] = -numerator / denominator;
        }
        
        // Calculate delta_p and perform collision detection and response
        double corr_q = Kernels::poly6(Constants::corr_q * Constants::h * dvec3{1.0, 0.0, 0.0}, Constants::h);
        for (int i = 0; i < scene.particles.size(); ++i) {
            delta_pos[i] = dvec3{0.0, 0.0, 0.0};
            for (int j = neighbor_offsets[i]; j < neighbor_offsets[i] + neighbor_sizes[i]; j++) {
                int neighbor_id = neighbors[j];
                double corr_kernel = Kernels::poly6(scene.particles[i].new_pos - scene.particles[neighbor_id].new_pos, Constants::h); 
                double corr = -Constants::corr_k * std::pow(corr_kernel / corr_q, Constants::corr_n);
                dvec3 grad_W = Kernels::gradSpiky(scene.particles[i].new_pos - scene.particles[neighbor_id].new_pos, Constants::h);
                delta_pos[i] += Constants::mass * (lambdas[i] + lambdas[neighbor_id] + corr) * grad_W;
            }
            delta_pos[i] *= (1.0 / Constants::mass) * (1.0 / Constants::rest_density);
        }
        
        // Separate for loops for parallelization later 
        for (int i = 0; i < scene.particles.size(); ++i) {
            // Update positions
            scene.particles[i].new_pos += delta_pos[i];
        }


        for (int i = 0; i < scene.particles.size(); ++i) {
            // Should we zero out velocities? 
            for (int j = 0; j < 3; j++) {
                if (scene.particles[i].new_pos[j] < bbox_mins[j] + Constants::radius) {
                    scene.particles[i].new_pos[j] = bbox_mins[j] + Constants::radius;
                } else if (scene.particles[i].new_pos[j] > bbox_maxs[j] - Constants::radius) {
                    scene.particles[i].new_pos[j] = bbox_maxs[j] - Constants::radius;
                }
            }
        }
    }

    for (auto& p : scene.particles) {
        // should we damp the velocities? 
        p.vel = Constants::damping * (p.new_pos - p.pos) / elapsed;
        p.pos = p.new_pos;
    }

    // TODO: vorticity 

    // XSPH viscosity
    for (int i = 0; i < scene.particles.size(); ++i) {
        densities[i] = compute_density(scene, i);
    }
    for (int i = 0; i < scene.particles.size(); ++i) {
        delta_vel[i] = dvec3{0.0, 0.0, 0.0};
        for (int j = neighbor_offsets[i]; j < neighbor_offsets[i] + neighbor_sizes[i]; j++) {
            int neighbor_id = neighbors[j];
            dvec3 vel = scene.particles[neighbor_id].vel - scene.particles[i].vel;
            double density = densities[neighbor_id];
            // delta_vel[i] += (Constants::xsph_c / density) * vel * Kernels::poly6(scene.particles[i].new_pos - scene.particles[neighbor_id].new_pos, Constants::h);
            delta_vel[i] += (Constants::mass / density) * vel * Kernels::poly6(scene.particles[i].new_pos - scene.particles[neighbor_id].new_pos, Constants::h);
        }
        delta_vel[i] *= Constants::xsph_c;
    }
    for (int i = 0; i < scene.particles.size(); ++i) {
        scene.particles[i].vel += delta_vel[i];
    }    
}
