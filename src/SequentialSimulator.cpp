#include "Simulator.hpp"
#include "Constants.hpp"

SequentialSimulator::SequentialSimulator() {}

void SequentialSimulator::update(double elapsed, Scene& scene) {
    // Apply gravity on all particles
    for (auto& p : scene.particles) {
        p.vel += elapsed * Constants::g;
        p.pos += elapsed * p.vel;
    }
}