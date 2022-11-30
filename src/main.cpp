#include "Constants.hpp"
#include "Scene.hpp"
#include "Simulator.hpp"
#include <iostream>
#include "Renderer.hpp"
#include <memory>
#include <GLFW/glfw3.h>
#include "Renderer.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    cout << "Starting main" << endl;
    // todo: parse arguments (parallel/sequential, scene file, benchmark/visualization, num iterations)
    bool parallel = false;
    string scene_file = "../scenes/one_particle.txt";
    bool benchmark = false;
    int num_iterations = 1000000;

    cout << "Loading scene" << endl;
    Scene scene(scene_file);
    cout << scene.particles.size() << " particles loaded" << endl;

    cout << "Initiating sim" << endl;
    unique_ptr<Simulator> sim;
    if (parallel) sim.reset(new ParallelSimulator(scene));
    else sim.reset(new SequentialSimulator(scene));


    cout << "Preparing graphics" << endl;
    GLFWwindow* window;
    if (!glfwInit()) {
        cout << "glfwInit() failed" << endl;
        return -1;
    }
    window = glfwCreateWindow(1280, 720, "Epic Fluid Simulation", NULL, NULL);
    if (!window) {
        cout << "Unable to create window" << endl;
        glfwTerminate();
        return -1;
    }

    cout << "Preparing renderer" << endl;
    Renderer renderer(window);

    cout << "Starting main loop" << endl;
    /* Loop until the user closes the window */
    for (int i = 0; i < num_iterations; i++) {
        if (glfwWindowShouldClose(window)) break;

        /* Render here */
        sim->update(Constants::dt, scene);
        renderer.draw(scene);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
