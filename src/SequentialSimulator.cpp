#include "Simulator.hpp"
#include "Constants.hpp"

SequentialSimulator::SequentialSimulator() {}

void SequentialSimulator::update(double elapsed, Scene& scene) {
    for (auto& p : scene.particles) {
        // Apply forces
        p.vel += elapsed * Constants::g;

        // Predict new position
        p.new_pos = p.pos + elapsed * p.vel;
        // Collisions with box
        for (int i = 0; i < 3; i++) {
            if (p.new_pos[i] < scene.bbox_mins[i] + Constants::radius) {
                p.vel[i] = 0;
                p.new_pos[i] = scene.bbox_mins[i] + Constants::radius;
            } else if (p.new_pos[i] > scene.bbox_maxs[i] - Constants::radius) {
                p.vel[i] = 0;
                p.new_pos[i] = scene.bbox_maxs[i] - Constants::radius;
            }
        }
    }

    for (auto& p : scene.particles) {
        p.vel = 1.0 / elapsed * (p.new_pos - p.pos);
        // TODO: apply vorticity confinement and XSPH viscosity
        p.pos = p.new_pos;
    }
}