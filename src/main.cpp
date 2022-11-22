#include "Constants.hpp"
#include "Scene.hpp"
#include "Simulator.hpp"
#include "Renderer.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    // todo: parse arguments (parallel/sequential, scene file, benchmark/visualization, num iterations)
    bool parallel = false;
    string scene_file = "one_particle.txt";
    bool benchmark = false;
    int num_iterations = 100;

    Simulator sim;
    if (parallel) sim = ParallelSimulator();
    else sim = SequentialSimulator();

    Scene scene(scene_file);

    /*Renderer renderer;

    for (int i = 0; i < num_iterations; i++) {
        sim.update(Constants::dt, scene);
        if (!benchmark) {
            renderer.draw(scene);
        }
    }*/
}
