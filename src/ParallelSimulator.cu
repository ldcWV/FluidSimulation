#include "Simulator.hpp"
#include "Particle.hpp"
#include "Constants.hpp"
#include "cuda_runtime.h"
#include "cuda.h"
#include <thrust/scan.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

__host__ ParallelSimulator::ParallelSimulator(const Scene& scene) {
    _n = scene.particles.size();
    cudaMalloc(&n, sizeof(size_t));
    cudaMemset(&n, _n, sizeof(size_t));   
    cudaMalloc(&threads, sizeof(size_t));
    cudaMemset(&threads, Constants::threads_per_block, sizeof(size_t));
    _blocks = (_n + Constants::threads_per_block - 1) / Constants::threads_per_block;
    cudaMalloc(&blocks, sizeof(size_t));
    cudaMemset(&blocks, _blocks, sizeof(size_t));

    cudaMalloc(&bbox_mins, sizeof(glm::dvec3));
    cudaMalloc(&bbox_maxs, sizeof(glm::dvec3));
    cudaMemcpy(bbox_mins, &scene.bbox_mins, sizeof(glm::dvec3), cudaMemcpyHostToDevice);
    cudaMemcpy(bbox_maxs, &scene.bbox_maxs, sizeof(glm::dvec3), cudaMemcpyHostToDevice);

    cudaMalloc(&grid_width, sizeof(size_t));
    cudaMalloc(&grid_height, sizeof(size_t));
    cudaMalloc(&grid_length, sizeof(size_t));
    size_t width = (scene.bbox_maxs.x - scene.bbox_mins.x) / Constants::h + 1;
    size_t height = (scene.bbox_maxs.y - scene.bbox_mins.y) / Constants::h + 1;
    size_t length = (scene.bbox_maxs.z - scene.bbox_mins.z) / Constants::h + 1;
    cudaMemset(grid_width, width, sizeof(size_t));
    cudaMemset(grid_height, height, sizeof(size_t));
    cudaMemset(grid_length, length, sizeof(size_t));

    cudaMalloc((void**)&lambdas, sizeof(double) * n);
    cudaMalloc((void**)&densities, sizeof(double) * n);
    cudaMalloc((void**)&particles, sizeof(Particle) * n);
    cudaMalloc((void**)&sorted_particles, sizeof(Particle) * n);
    cudaMalloc((void**)&neighbor_starts, sizeof(int) * n);

    _total_cells = width * height * length; // overflow?
    cudaMalloc(&total_cells, sizeof(size_t));
    cudaMemset(&total_cells, _total_cells, sizeof(size_t)); 
    cudaMalloc((void**)&bins, sizeof(int) * total_cells);
    cudaMalloc((void**)&prefix_bins, sizeof(int) * total_cells);

    cudaMalloc((void**)&delta_pos, sizeof(glm::dvec3) * n);
    cudaMalloc((void**)&delta_vel, sizeof(glm::dvec3) * n);
}

__host__ ParallelSimulator::~ParallelSimulator() {
    cudaFree(n);
    cudaFree(threads);
    cudaFree(blocks);
    cudaFree(total_cells);
    cudaFree(bbox_mins);
    cudaFree(bbox_maxs);
    cudaFree(grid_width);
    cudaFree(grid_height);
    cudaFree(grid_length);
    cudaFree(lambdas);
    cudaFree(densities);
    cudaFree(particles);
    cudaFree(sorted_particles);
    cudaFree(bins);
    cudaFree(prefix_bins);
    cudaFree(delta_pos);
    cudaFree(delta_vel);
}

__host__ void ParallelSimulator::reset() {
    // Consider using kernels to zero out memory 
    cudaMemset(lambdas, 0, sizeof(double) * _n);
    cudaMemset(densities, 0, sizeof(double) * _n);
    cudaMemset(neighbor_starts, 0, sizeof(int) * _n)
    cudaMemset(bins, 0, sizeof(int) * _total_cells);
    cudaMemset(prefix_bins, 0, sizeof(int) * _total_cells);
    cudaMemset(delta_pos, 0, sizeof(glm::dvec3) * _n);
    cudaMemset(delta_vel, 0, sizeof(glm::dvec3) * _n);
}

__device__ double poly6(const glm::dvec3& r, const double h) {
    double r_mag = glm::length(r);
    if (r_mag > h) return 0.0;
    return 315.0 / (64 * Constants::pi * pow(h, 9)) * pow(h * h - r_mag * r_mag, 3);
}

__device__ glm::dvec3 grad_spiky(const glm::dvec3& r, const double h) {
    double r_mag = glm::length(r);
    if (r_mag > h) return glm::dvec3{0, 0, 0};
    return -45 / (Constants::pi * pow(h, 6) * std::max(r_mag, 1e-24)) * pow(h - r_mag, 2) * r;
}

__device__ glm::ivec3 get_cell_coords(glm::dvec3 pos) {
    ivec3 res {
        (pos.x - bbox_mins.x) / Constants::h,
        (pos.y - bbox_mins.y) / Constants::h,
        (pos.z - bbox_mins.z) / Constants::h
    };
    return clamp(res, ivec3{0, 0, 0}, ivec3{grid_width-1, grid_height-1, grid_length-1});
}

__device__ int get_cell_idx(glm::ivec3 coords) {
    return coords.x * grid_height * grid_length + coords.y * grid_length + coords.z;
    // return coords.x + coords.y * grid_width + coords.z * grid_width * grid_height;
} 

__global__ void compute_bins_kernel(Particle *particles, int *bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    glm::ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    atomicAdd(&bins[cell_idx], 1);
}

__host__ void ParallelSimulator::compute_bins() {
    compute_bins_kernel<<<blocks, threads>>>(particles, bins, n);
    cudaDeviceSynchronize();
}

__host__ void ParallelSimulator::compute_prefix_bins() {
    thrust::exclusive_scan(thrust::device, bins, bins + total_cells, prefix_bins);
    cudaDeviceSynchronize();
}

__global__ void compute_sorted_particles_kernel(Particle *particles, int *prefix_bins, Particle *sorted_particles, n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    glm::ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    // Returns the value of prefix_bins[cell_idx] before the addition
    // THIS MIGHT BE BAD if atomic Add is not actually atomic 
    // cause multiple threads might get the same value of prefix_bins[cell_idx]
    int sorted_idx = atomicAdd(&prefix_bins[cell_idx], 1);
    sorted_particles[sorted_idx] = p;
}

__host__ void ParallelSimulator::compute_sorted_particles() {
    compute_sorted_particles_kernel<<<blocks, threads>>>(particles, prefix_bins, sorted_particles, n);
    cudaDeviceSynchronize();
}

__global__ void compute_neighbor_starts_kernel(Particle *sorted_particles, int *neighbor_starts, int *prefix_bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    glm::ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    if (cell_idx == 0) {
        neighbor_starts[idx] = 0;
    } else {
        neighbor_starts[idx] = prefix_bins[cell_idx - 1];
    }
}

__host__ void ParallelSimulator::compute_neighbor_starts() {
    compute_neighbor_starts_kernel<<<blocks, threads>>>(sorted_particles, neighbor_starts, prefix_bins, n);
    cudaDeviceSynchronize();
}

__host__ void ParallelSimulator::recompute_neighbors() {
    compute_bins();
    compute_prefix_bins();
    compute_sorted_particles();
    compute_neighbor_starts();
}

__device__ void compute_constraint(double *densities, int particle_id) {
    return densities[particle_id] / Constants::rest_density - 1.0;
}

__device__ glm::dvec3 compute_grad_constraint(Particle* sorted_particles, int neighbor_start, int neighbor_end, int constraint_id, int grad_id) {
    const Particle &constraint_particle = sorted_particles[constraint_id];
    const glm::dvec3& constraint_pos = constraint_particle.new_pos;
    if (constraint_id == grad_id) {
        glm::dvec3 res{0.0, 0.0, 0.0};
        for (int ni = neighbor_start; ni < neighbor_end; ni++) {
            const Particle &neighbor = sorted_particles[ni];
            const glm::dvec3& neighbor_pos = neighbor.new_pos;
            res += Constants::mass * grad_spiky(constraint_pos - neighbor_pos, Constants::h);
        }
        return res / Constants::rest_density;
    } else {
        return -Constants::mass * grad_spiky(constraint_pos - sorted_particles[grad_id].new_pos, Constants::h) / Constants::rest_density;
    }
}

__global__ void compute_lambdas_kernel(Particle *sorted_particles, double *densities, double *lambdas, int *neighbor_starts, int *bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    glm::ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    double numerator = compute_constraint(densities, idx);
    double denominator = 0.0;
    int neighbor_start = neighbor_starts[idx];
    int neighbor_end = neighbor_start + bins[cell_idx];
    for (int ni = neighbor_start; ni < neighbor_end; ni++) {
        glm::dvec3 grad = compute_grad_constraint(sorted_particles, neighbor_start, neighbor_end, idx, ni);
        denominator += glm::dot(grad, grad) / Constants::mass;
    }
    denominator += Constants::epsilon;
    lambdas[idx] = -numerator / denominator;
}

__host__ void ParallelSimulator::compute_lambdas() {
    compute_lambdas_kernel<<<blocks, threads>>>(sorted_particles, densities, lambdas, neighbor_starts, prefix_bins, n);
    cudaDeviceSynchronize();
}

__global__ void compute_delta_positions_kernel(Particle *sorted_particles, double *delta_pos, double *lambdas, int *neighbor_starts, int *bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    glm::ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);

    double corr_q = poly6(Constants::corr_q * Constants::h * dvec3{1.0, 1.0, 1.0}, Constants::h);
    delta_pos[idx] = glm::dvec3{0.0, 0.0, 0.0};

    int neighbor_start = neighbor_starts[idx];
    int neighbor_end = neighbor_start + bins[cell_idx];
    for (int ni = neighbor_start; ni < neighbor_end; ni++) {
        const Particle &neighbor = sorted_particles[ni];
        double corr_kernel = poly6(p.new_pos - neighbor.new_pos, Constants::h); 
        double corr = -Constants::corr_k * std::pow(corr_kernel / corr_q, Constants::corr_n);
        glm::dvec3 grad_W = grad_spiky(p.new_pos - neighbor.new_pos, Constants::h);
        glm::dvec3 grad_W2 = compute_grad_constraint(sorted_particles, neighbor_start, neighbor_end, idx, ni);
        delta_pos[idx] += Constants::mass * (lambdas[idx] + lambdas[ni] + corr) * grad_W;
    }
    delta_pos[idx] *= (1.0 / Constants::mass) * (1.0 / Constants::rest_density);
}

__host__ void ParallelSimulator::compute_delta_positions() {
    compute_delta_positions<<<blocks, threads>>>(sorted_particles, delta_pos, lambdas, neighbor_starts, prefix_bins, n)
    cudaDeviceSynchronize();
}


__global__ void update_positions_kernel(Particle *sorted_particles, glm::dvec3 *delta_pos, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    p.new_pos += delta_pos[idx];
}

__host__ void ParallelSimulator::update_positions() {
    update_positions_kernel(sorted_particles, delta_pos, n);
    cudaDeviceSynchronize();
}

__global__ void update_collisions_kernel(Particle *sorted_particles, glm::dvec3 bbox_mins, glm::dvec3 bbox_maxs, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    for (int j = 0; j < 3; j++) {
        if (p.new_pos[j] < bbox_mins[j] + Constants::radius) {
            p.new_pos[j] = bbox_mins[j] + Constants::radius;
        } else if (p.new_pos[j] > bbox_maxs[j] - Constants::radius) {
            p.new_pos[j] = bbox_maxs[j] - Constants::radius;
        }
    }
}

__host__ void ParallelSimulator::update_collisions() {
    update_collisions_kernel(sorted_particles, bbox_mins, bbox_maxs, n);
    cudaDeviceSynchronize();
}

__host__ void ParallelSimulator::simulate() {
    compute_densities(); // to use in lambdas and delta positions
    compute_lambdas();
    compute_delta_positions();
    update_positions();
    update_collisions();
}

__global__ void compute_velocities_kernel(Particle *particles, double elapsed, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    p.vel += elapsed * Constants::g;
    p.new_pos = p.pos + elapsed * p.vel;
}

__host__ void ParallelSimulator::compute_velocities(double elapsed) {
    compute_velocities_kernel<<<blocks, threads>>>(particles, elapsed, n);
    cudaDeviceSynchronize();
}

__global__ void compute_velocities_and_positions_kernel(Particle *sorted_particles, double elapsed, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    p.vel = Constants::damping * (p.new_pos - p.pos) / elapsed;
    p.pos = p.new_pos;
}

__host__ void ParallelSimulator::compute_velocities_and_positions(double elapsed) {
    compute_velocities_and_positions_kernel<<<blocks, threads>>>(sorted_particles, elapsed, n);
    cudaDeviceSynchronize();
}

__global__ void compute_densities_kernel(Particle *sorted_particles, double *densities, int *neighbor_starts, int *bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    glm::ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    double density = 0.0;
    for (int ni = neighbor_starts[idx]; ni < neighbor_starts[idx] + bins[cell_idx]; ni++) {
        Particle &neighbor = sorted_particles[ni];
        density += Constant::mass * poly6(p.new_pos - neighbor.new_pos, Constants::h);
    }
    densities[idx] = density / Constants::rest_density - 1.0;
}

__host__ void ParallelSimulator::compute_densities() {
    compute_densities_kernel<<<blocks, threads>>>(sorted_particles, densities, neighbor_starts, bins, n);
    cudaDeviceSynchronize();
}

__global__ void xsph_viscosity_kernel(Particle *sorted_particles, glm::dvec3 *delta_vel, double *densities, int *neighbor_starts, int *bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    glm::dvec3 &delta_vel = delta_vel[idx];
    glm::ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    double delta_vel = 0.0;
    for (int ni = neighbor_starts[idx]; ni < neighbor_starts[idx] + bins[cell_idx]; ni++) {
        Particle &neighbor = sorted_particles[ni];
        glm::dvec3 vel = neighbor.vel - p.vel;
        double density = densities[ni];
        delta_vel += (Constants::mass / density) * poly6(p.new_pos - neighbor.new_pos, Constants::h);
    }
    delta_vel[idx] = delta_vel;
}

__host__ void ParallelSimulator::xsph_viscosity() {
    xsph_viscosity_kernel<<<blocks, threads>>>(sorted_particles, delta_vel, densities, neighbor_starts, bins, n);
    cudaDeviceSynchronize();
}

__global__ void update_velocities_kernel(Particle *sorted_particles, glm::dvec3 *delta_vel, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = sorted_particles[idx];
    p.vel += delta_vel[idx];
}

__host__ void ParallelSimulator::update_velocities() {
    update_velocities_kernel<<<blocks, threads>>>(sorted_particles, delta_vel, n);
    cudaDeviceSynchronize();
}

__host__ void ParallelSimulator::update(double elapsed, Scene& scene) {
    reset();

    // Copy scene particles from host to device memory 
    cudaMemcpy(particles, scene.particles.data(), sizeof(Particle) * _n, cudaMemcpyHostToDevice);

    // Initial forces 
    compute_velocities(elapsed);

    // Recompute neighbors 
    recompute_neighbors();

    // Simulation 
    for (int iter = 0; iter < Constants::solver_iterations; iter++) {
        simulate(); 
    }

    // Post processing (update velocities, positions, XSPH, vorticity)
    // TODO: vorticity
    compute_densities();
    compute_velocities_and_positions(elapsed);
    xsph_visocisty();
    update_velocities();

    // Copy particles back to host 
    cudaMemcpy(scene.particles.data(), sorted_particles, sizeof(Particle) * _n, cudaMemcpyDeviceToHost);
    return;
}