#include "Simulator.hpp"
#include "Particle.hpp"
#include "Constants.hpp"
#include "cuda_runtime.h"
#include "cuda.h"
#include <thrust/scan.h>
#include <GLFW/glfw3.h>
// #include <thrust/device_ptr.h>
// #include <thrust/device_malloc.h>
// #include <thrust/device_free.h>

#define cudaCheckError(ans) cudaAssert((ans), __FILE__, __LINE__);
inline void cudaAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {   
        fprintf(stderr, "CUDA Error: %s at %s:%d\n",
        cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

using namespace glm;

// consider moving Constants to here
struct GlobalConstants {
    dvec3 bbox_mins;
    dvec3 bbox_maxs; 
    size_t grid_width;
    size_t grid_height;
    size_t grid_length; 
    dvec3 g;
    double eps;
    double mass;
    double h;
    double radius;
    int solver_iterations;
    double pi;
    double rest_density;
    double corr_q;
    double corr_k;
    int corr_n;
    float xsph_c;
    double damping;
    int threads_per_block;
};

__constant__ GlobalConstants GC; 

__host__ ParallelSimulator::ParallelSimulator(const Scene& scene) {
    _n = scene.particles.size();
    cudaMalloc((void**)&n, sizeof(size_t));
    cudaMemcpy((void**)&n, &_n, sizeof(size_t), cudaMemcpyHostToDevice);   
    _blocks = (_n + Constants::threads_per_block - 1) / Constants::threads_per_block;
    _threads = Constants::threads_per_block;

    cudaMalloc((void**)&bbox_mins, sizeof(dvec3));
    cudaMalloc((void**)&bbox_maxs, sizeof(dvec3));
    cudaMemcpy(&bbox_mins, &scene.bbox_mins, sizeof(dvec3), cudaMemcpyHostToDevice);
    cudaMemcpy(&bbox_maxs, &scene.bbox_maxs, sizeof(dvec3), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&grid_width, sizeof(size_t));
    cudaMalloc((void**)&grid_height, sizeof(size_t));
    cudaMalloc((void**)&grid_length, sizeof(size_t));
    size_t width = (scene.bbox_maxs.x - scene.bbox_mins.x) / Constants::h + 1;
    size_t height = (scene.bbox_maxs.y - scene.bbox_mins.y) / Constants::h + 1;
    size_t length = (scene.bbox_maxs.z - scene.bbox_mins.z) / Constants::h + 1;
    cudaMemcpy(&grid_width, &width, sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(&grid_height, &height, sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(&grid_length, &length, sizeof(size_t), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&lambdas, sizeof(double) * _n);
    cudaMalloc((void**)&densities, sizeof(double) * _n);
    cudaMalloc((void**)&particles, sizeof(Particle) * _n);
    cudaMalloc((void**)&neighbors, sizeof(int) * 100000000);
    cudaMalloc((void**)&neighbor_starts, sizeof(int) * _n);
    cudaMalloc((void**)&neighbor_sizes, sizeof(int) * _n);

    _total_cells = width * height * length; // overflow?
    cudaMalloc((void**)&total_cells, sizeof(size_t));
    cudaMemcpy(&total_cells, &_total_cells, sizeof(size_t), cudaMemcpyHostToDevice); 
    cudaMalloc((void**)&bins, sizeof(int) * _total_cells);
    cudaMalloc((void**)&prefix_bins, sizeof(int) * _total_cells);
    cudaMalloc((void**)&grid_starts, sizeof(int) * _total_cells);
    cudaMalloc((void**)&grid, sizeof(int) * _n);

    cudaMalloc((void**)&delta_pos, sizeof(dvec3) * _n);
    cudaMalloc((void**)&delta_vel, sizeof(dvec3) * _n);

    GlobalConstants _GC;
    _GC.bbox_mins = scene.bbox_mins;
    _GC.bbox_maxs = scene.bbox_maxs;
    _GC.grid_width = width;
    _GC.grid_height = height;
    _GC.grid_length = length;
    _GC.g = Constants::g;
    _GC.eps = Constants::eps;
    _GC.mass = Constants::mass;
    _GC.h = Constants::h;
    _GC.radius = Constants::radius;
    _GC.solver_iterations = Constants::solver_iterations;
    _GC.pi = Constants::pi;
    _GC.rest_density = Constants::rest_density;
    _GC.corr_q = Constants::corr_q;
    _GC.corr_k = Constants::corr_k;
    _GC.corr_n = Constants::corr_n;
    _GC.xsph_c = Constants::xsph_c;
    _GC.damping = Constants::damping;
    _GC.threads_per_block = Constants::threads_per_block;

    cudaMemcpyToSymbol(GC, &_GC, sizeof(GlobalConstants));
}

__host__ ParallelSimulator::~ParallelSimulator() {
    cudaFree(&n);
    cudaFree(&total_cells);
    cudaFree(&bbox_mins);
    cudaFree(&bbox_maxs);
    cudaFree(&grid_width);
    cudaFree(&grid_height);
    cudaFree(&grid_length);
    cudaFree(lambdas);
    cudaFree(densities);
    cudaFree(particles);
    cudaFree(neighbors);
    cudaFree(bins);
    cudaFree(prefix_bins);
    cudaFree(delta_pos);
    cudaFree(delta_vel);
    cudaFree(neighbor_sizes);
    cudaFree(neighbor_starts);
    cudaFree(grid);
}

__host__ void ParallelSimulator::reset() {
    // Consider using kernels to zero out memory 
    cudaMemset(lambdas, 0, sizeof(double) * _n);
    cudaMemset(densities, 0, sizeof(double) * _n);
    cudaMemset(neighbor_sizes, 0, sizeof(int) * _n);
    cudaMemset(bins, 0, sizeof(int) * _total_cells);
    cudaMemset(prefix_bins, 0, sizeof(int) * _total_cells);
    cudaMemset(delta_pos, 0, sizeof(dvec3) * _n);
    cudaMemset(delta_vel, 0, sizeof(dvec3) * _n);
}

__device__ double poly6(const dvec3& r, const double h) {
    double r_mag = length(r);
    if (r_mag > h) return 0.0;
    double powh3 = h*h*h;
    double powh9 = powh3*powh3*powh3;
    double h_rmag = h*h - r_mag*r_mag;
    double h_rmag3 = h_rmag*h_rmag*h_rmag;
    return 315.0 / (64 * GC.pi * powh9) * h_rmag3;
}

__device__ dvec3 grad_spiky(const dvec3& r, const double h) {
    double r_mag = length(r);
    if (r_mag > h) return dvec3{0, 0, 0};
    double powh2 = h*h;
    double powh6 = powh2*powh2*powh2;
    double dist2 = dot(h-r_mag, h-r_mag);
    return -45 / (GC.pi * powh6 * max(r_mag, 1e-24)) * dist2 * r;
}

__device__ ivec3 get_cell_coords(dvec3 pos) {
    ivec3 res {
        (pos.x - GC.bbox_mins.x) / GC.h,
        (pos.y - GC.bbox_mins.y) / GC.h,
        (pos.z - GC.bbox_mins.z) / GC.h
    };
    return clamp(res, ivec3{0, 0, 0}, ivec3{GC.grid_width-1, GC.grid_height-1, GC.grid_length-1});
}

__device__ int get_cell_idx(ivec3 coords) {
    return coords.x * GC.grid_height * GC.grid_length + coords.y * GC.grid_length + coords.z;
} 

__global__ void compute_bins_kernel(Particle *particles, int *bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    atomicAdd(&bins[cell_idx], 1);
}

__global__ void compute_sorted_grid_kernel(Particle *particles, int *prefix_bins, int *grid, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    int sorted_idx = atomicAdd(&prefix_bins[cell_idx], 1);
    grid[sorted_idx] = p.id;
}

__global__ void compute_grid_starts_kernel(Particle *particles, int *grid_starts, int *prefix_bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    if (cell_idx == 0) {
        grid_starts[cell_idx] = 0;
    } else {
        grid_starts[cell_idx] = prefix_bins[cell_idx - 1];
    }
}


__host__ void ParallelSimulator::recompute_grid() {
    compute_bins_kernel<<<_blocks, _threads>>>(particles, bins, _n);
    // cudaCheckError(cudaDeviceSynchronize());
    thrust::exclusive_scan(thrust::device, bins, bins + _total_cells, prefix_bins);
    // cudaCheckError(cudaDeviceSynchronize());
    compute_sorted_grid_kernel<<<_blocks, _threads>>>(particles, prefix_bins, grid, _n);
    // cudaCheckError(cudaDeviceSynchronize());
    compute_grid_starts_kernel<<<_blocks, _threads>>>(particles, grid_starts, prefix_bins, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__global__ void compute_neighbor_sizes_kernel(Particle *particles, int *neighbor_sizes, int *bins, int _n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= _n) return;
    Particle &p = particles[idx];
    ivec3 our_coords = get_cell_coords(p.new_pos);

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                ivec3 coords = our_coords + ivec3{dx, dy, dz};
                if (coords.x < 0 || coords.x >= GC.grid_width ||
                    coords.y < 0 || coords.y >= GC.grid_height ||
                    coords.z < 0 || coords.z >= GC.grid_length) continue;
                
                int cell_idx = get_cell_idx(coords);
                neighbor_sizes[idx] += bins[cell_idx];
            }
        }
    }
}


__global__ void compute_neighbors_kernel(Particle *particles, int *neighbors, int *neighbor_sizes, int *neighbor_starts, int *grid, int *grid_starts, int *bins, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    ivec3 our_coords = get_cell_coords(p.new_pos);
    // Hacky way to make ourself always the first neighbor
    neighbors[neighbor_starts[idx]] = idx;
    int current_neighbor = 1;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                ivec3 coords = our_coords + ivec3{dx, dy, dz};
                if (coords.x < 0 || coords.x >= GC.grid_width ||
                    coords.y < 0 || coords.y >= GC.grid_height ||
                    coords.z < 0 || coords.z >= GC.grid_length) continue;
                
                int cell_idx = get_cell_idx(coords);
                for (int j = grid_starts[cell_idx]; j < grid_starts[cell_idx] + bins[cell_idx]; ++j) {
                    if (grid[j] != idx) {
                        neighbors[neighbor_starts[idx] + current_neighbor] = grid[j];
                        current_neighbor++;
                    }
                }
            }
        }
    }
}

__host__ void ParallelSimulator::recompute_neighbors() {
    compute_neighbor_sizes_kernel<<<_blocks, _threads>>>(particles, neighbor_sizes, bins, _n);
    // cudaCheckError(cudaDeviceSynchronize());
    thrust::exclusive_scan(thrust::device, neighbor_sizes, neighbor_sizes + _n, neighbor_starts);
    // cudaCheckError(cudaDeviceSynchronize());
    compute_neighbors_kernel<<<_blocks, _threads>>>(particles, neighbors, neighbor_sizes, neighbor_starts, grid, grid_starts, bins, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__device__ double compute_constraint(double *densities, int particle_id) {
    return densities[particle_id] / GC.rest_density - 1.0;
}

__device__ dvec3 compute_grad_constraint(Particle* particles, int constraint_id, int grad_id) {
    const Particle &constraint_particle = particles[constraint_id];
    const dvec3& constraint_pos = constraint_particle.new_pos;
    return -GC.mass * grad_spiky(constraint_pos - particles[grad_id].new_pos, GC.h) / GC.rest_density;
}

__device__ dvec3 compute_grad_constraint_self(Particle* particles, int* neighbors, int neighbor_start, int neighbor_end, int constraint_id) {
    const Particle &constraint_particle = particles[constraint_id];
    const dvec3& constraint_pos = constraint_particle.new_pos;
    dvec3 res{0.0, 0.0, 0.0};
    for (int ni = neighbor_start; ni < neighbor_end; ni++) {
        int neighbor_id = neighbors[ni];
        Particle& neighbor = particles[neighbor_id];
        const dvec3& neighbor_pos = neighbor.new_pos;
        res += GC.mass * grad_spiky(constraint_pos - neighbor_pos, GC.h);
    }
    return res / GC.rest_density;
}

__global__ void compute_lambdas_kernel(Particle *particles, int *neighbors, double *densities, double *lambdas, int *neighbor_starts, int *neighbor_sizes, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    int neighbor_start = neighbor_starts[idx];
    int neighbor_end = neighbor_start + neighbor_sizes[idx];
    double numerator = compute_constraint(densities, idx);
    dvec3 grad_self  = compute_grad_constraint_self(particles, neighbors, neighbor_start, neighbor_end, idx);
    double denominator = dot(grad_self, grad_self) / GC.mass;
    for (int ni = neighbor_start + 1; ni < neighbor_end; ni++) {
        int neighbor_id = neighbors[ni];
        Particle& neighbor = particles[neighbor_id];
        dvec3 grad = compute_grad_constraint(particles, idx, neighbor.id);
        denominator += dot(grad, grad) / GC.mass;
    }
    denominator += GC.eps;
    lambdas[idx] = -numerator / denominator;
}

__host__ void ParallelSimulator::compute_lambdas() {
    compute_lambdas_kernel<<<_blocks, _threads>>>(particles, neighbors, densities, lambdas, neighbor_starts, neighbor_sizes, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__global__ void compute_delta_positions_kernel(Particle *particles, int *neighbors, dvec3 *delta_pos, double *lambdas, int *neighbor_starts, int *neighbor_sizes, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);

    double corr_q = poly6(GC.corr_q * GC.h * dvec3{1.0, 0.0, 0.0}, GC.h);
    delta_pos[idx] = dvec3{0.0, 0.0, 0.0};

    int neighbor_start = neighbor_starts[idx];
    int neighbor_end = neighbor_start + neighbor_sizes[idx];
    for (int ni = neighbor_start; ni < neighbor_end; ni++) {
        int neighbor_id = neighbors[ni];
        Particle &neighbor = particles[neighbor_id];
        double corr_kernel = poly6(p.new_pos - neighbor.new_pos, GC.h); 
        double ratio = corr_kernel / corr_q;
        double ratio2 = ratio*ratio;
        double ratio4 = ratio2*ratio2;
        double corr = -GC.corr_k * ratio4;
        dvec3 grad_W = grad_spiky(p.new_pos - neighbor.new_pos, GC.h);
        delta_pos[idx] += GC.mass * (lambdas[idx] + lambdas[neighbor.id] + corr) * grad_W;
    }
    delta_pos[idx] *= (1.0 / GC.mass) * (1.0 / GC.rest_density);
}


__host__ void ParallelSimulator::compute_delta_positions() {
    compute_delta_positions_kernel<<<_blocks, _threads>>>(particles, neighbors, delta_pos, lambdas, neighbor_starts, neighbor_sizes, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}


__global__ void update_positions_kernel(Particle *particles, dvec3 *delta_pos, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    p.new_pos += delta_pos[idx];
}

__host__ void ParallelSimulator::update_positions() {
    update_positions_kernel<<<_blocks, _threads>>>(particles, delta_pos, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__global__ void update_collisions_kernel(Particle *particles, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    for (int j = 0; j < 3; j++) {
        if (p.new_pos[j] < GC.bbox_mins[j] + GC.radius) {
            p.new_pos[j] = GC.bbox_mins[j] + GC.radius;
        } else if (p.new_pos[j] > GC.bbox_maxs[j] - GC.radius) {
            p.new_pos[j] = GC.bbox_maxs[j] - GC.radius;
        }
    }
}

__host__ void ParallelSimulator::update_collisions() {
    update_collisions_kernel<<<_blocks, _threads>>>(particles, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__host__ void ParallelSimulator::simulate() {
    double t0 = glfwGetTime();

    compute_densities(); // to use in lambdas and delta positions
    cudaCheckError(cudaDeviceSynchronize());
    double t1 = glfwGetTime();

    compute_lambdas();
    cudaCheckError(cudaDeviceSynchronize());
    double t2 = glfwGetTime();

    compute_delta_positions();
    cudaCheckError(cudaDeviceSynchronize());
    double t3 = glfwGetTime();

    update_positions();
    cudaCheckError(cudaDeviceSynchronize());
    double t4 = glfwGetTime();

    update_collisions();
    cudaCheckError(cudaDeviceSynchronize());
    double t5 = glfwGetTime();

    std::cout << t1 - t0 << " " << t2-t1 << " " << t3-t2 << " " << t4-t3 << " " << t5-t4 << std::endl;
}

__global__ void compute_velocities_kernel(Particle *particles, double elapsed, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    p.vel += elapsed * GC.g;
    p.new_pos = p.pos + elapsed * p.vel;
}

__host__ void ParallelSimulator::compute_velocities(double elapsed) {
    compute_velocities_kernel<<<_blocks, _threads>>>(particles, elapsed, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__global__ void compute_velocities_and_positions_kernel(Particle *particles, double elapsed, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    p.vel = GC.damping * (p.new_pos - p.pos) / elapsed;
    p.pos = p.new_pos;
}

__host__ void ParallelSimulator::compute_velocities_and_positions(double elapsed) {
    compute_velocities_and_positions_kernel<<<_blocks, _threads>>>(particles, elapsed, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__global__ void compute_densities_kernel(Particle *particles, int *neighbors, double *densities, int *neighbor_starts, int *neighbor_sizes, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    double density = 0.0;
    for (int ni = neighbor_starts[idx]; ni < neighbor_starts[idx] + neighbor_sizes[idx]; ni++) {
        int neighbor_id = neighbors[ni];
        Particle &neighbor = particles[neighbor_id];
        density += GC.mass * poly6(p.new_pos - neighbor.new_pos, GC.h);
    }
    densities[idx] = density;
}

__host__ void ParallelSimulator::compute_densities() {
    compute_densities_kernel<<<_blocks, _threads>>>(particles, neighbors, densities, neighbor_starts, neighbor_sizes, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__global__ void xsph_viscosity_kernel(Particle *particles, int *neighbors, dvec3 *delta_vel, double *densities, int *neighbor_starts, int *neighbor_sizes, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    ivec3 coords = get_cell_coords(p.new_pos);
    int cell_idx = get_cell_idx(coords);
    delta_vel[idx] = dvec3{0.0, 0.0, 0.0};
    for (int ni = neighbor_starts[idx]; ni < neighbor_starts[idx] + neighbor_sizes[idx]; ni++) {
        int neighbor_id = neighbors[ni];
        Particle &neighbor = particles[neighbor_id];
        dvec3 vel = neighbor.vel - p.vel;
        double density = densities[neighbor.id];
        delta_vel[idx] += (GC.mass / density) * vel * poly6(p.new_pos - neighbor.new_pos, GC.h);
    }
}

__host__ void ParallelSimulator::xsph_viscosity() {
    xsph_viscosity_kernel<<<_blocks, _threads>>>(particles, neighbors, delta_vel, densities, neighbor_starts, neighbor_sizes, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__global__ void update_velocities_kernel(Particle *particles, dvec3 *delta_vel, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    Particle &p = particles[idx];
    p.vel += delta_vel[idx];
}

__host__ void ParallelSimulator::update_velocities() {
    update_velocities_kernel<<<_blocks, _threads>>>(particles, delta_vel, _n);
    // cudaCheckError(cudaDeviceSynchronize());
}

__host__ void ParallelSimulator::update(double elapsed, Scene& scene) {
    double t0 = glfwGetTime();
    reset();
    cudaCheckError(cudaDeviceSynchronize());

    // Copy scene particles from host to device memory 
    double t1 = glfwGetTime();
    cudaMemcpy(particles, scene.particles.data(), sizeof(Particle) * _n, cudaMemcpyHostToDevice);
    cudaCheckError(cudaDeviceSynchronize());

    // Initial forces 
    double t2 = glfwGetTime();
    compute_velocities(elapsed);
    cudaCheckError(cudaDeviceSynchronize());

    // Recompute neighbors 
    double t3 = glfwGetTime();
    recompute_grid();
    recompute_neighbors();
    cudaCheckError(cudaDeviceSynchronize());

    // Simulation 
    double t4 = glfwGetTime();
    for (int iter = 0; iter < Constants::solver_iterations; iter++) {
        simulate(); 
    }
    cudaCheckError(cudaDeviceSynchronize());

    // Post processing (update velocities, positions, XSPH, vorticity)
    // TODO: vorticity
    double t5 = glfwGetTime();
    compute_densities();
    cudaCheckError(cudaDeviceSynchronize());
    double t6 = glfwGetTime();
    compute_velocities_and_positions(elapsed);
    cudaCheckError(cudaDeviceSynchronize());
    double t7 = glfwGetTime();
    xsph_viscosity();
    cudaCheckError(cudaDeviceSynchronize());
    double t8 = glfwGetTime();
    update_velocities();
    cudaCheckError(cudaDeviceSynchronize());

    // Copy particles back to host 
    double t9 = glfwGetTime();
    cudaMemcpy(scene.particles.data(), particles, sizeof(Particle) * _n, cudaMemcpyDeviceToHost);
    cudaCheckError(cudaDeviceSynchronize());
    double t10 = glfwGetTime();

    std::cout << "0 - 4: " << t1-t0 << " " << t2-t1 << " " << t3-t2 << " " << t4-t3 << " " << t5-t4 << std::endl;
    std::cout << "5 - 10: " << t6-t5 << " " << t7-t6 << " " << t8-t7 << " " << t9-t8 << " " << t10-t9 << std::endl;
    return;
}