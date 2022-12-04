#include "Simulator.hpp"
#include "cuda_runtime.h"
#include "cuda.h"
#include <thrust/scan.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

ParallelSimulator::ParallelSimulator(const Scene& scene) {


}

ParallelSimulator::~ParallelSimulator() {

}

void ParallelSimulator::update(double elapsed, Scene& scene) {
    return;
}

glm::ivec3 ParallelSimulator::get_cell_coords(glm::dvec3 pos) {
    return glm::ivec3{0, 0, 0};
}

__device__ __host__ int ParallelSimulator::get_cell_idx(glm::ivec3 coords) {
    return 0;
}
