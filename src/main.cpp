#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include "Constants.hpp"
#include "Scene.hpp"
#include "Simulator.hpp"
#include <iostream>
#include "Renderer.hpp"
#include <memory>
#include "Renderer.hpp"
#include <direct.h>
#include <thread>
#include <chrono>
#include <mutex>

using namespace std;

Scene renderer_scene;
mutex renderer_scene_mutex;
atomic_bool renderer_scene_updated = true;
atomic_bool renderer_done = false;

void render() {
    Renderer renderer;
    static Scene current_scene;

    while (true) {
        // atomically copy renderer_scene
        if (renderer_scene_updated) {
            renderer_scene_mutex.lock();
            renderer_scene_updated = false;
            current_scene = renderer_scene;
            renderer_scene_mutex.unlock();
        }

        // render the copied scene
        if (!renderer.draw(current_scene)) break;
    }

    renderer_done = true;
}

int main(int argc, char* argv[]) {
    cout << "Starting main" << endl;
    // todo: parse these arguments
    bool parallel = false;
    string scene_name = "10000_random_yuki";
    bool benchmark = false;
    int num_iterations = 1000000;
    bool save_replay = false;
    bool play_replay = true;

    if (save_replay && play_replay) {
        cout << "Cannot both save and play replay" << endl;
        return -1;
    }

    string scene_file = string(SCENE_DIR) + scene_name + ".txt";
    cout << "Loading scene from " << scene_file << endl;
    Scene scene(scene_file, false);
    cout << scene.particles.size() << " particles loaded" << endl;

    cout << "Initiating sim" << endl;
    unique_ptr<Simulator> sim;
    if (parallel) sim.reset(new ParallelSimulator(scene));
    else sim.reset(new SequentialSimulator(scene));

    string replay_folder = string(REPLAY_DIR) + scene_name;
    if (save_replay) {
        cout << "Creating replay folder at " << replay_folder << endl;
        _mkdir(replay_folder.c_str());
    }
    int replay_idx = 0;

    cout << "Starting renderer" << endl;
    renderer_scene = scene;
    thread render_thread(render);

    cout << "Starting main loop" << endl;
    for (int i = 0; i < num_iterations; i++) {
        string replay_file = string(replay_folder) + "/" + string(8 - std::min(8U, to_string(replay_idx).length()), '0') + to_string(replay_idx) + ".bin";
        replay_idx++;

        /* Render here */
        if (play_replay) {
            if (scene.load(replay_file, true) == -1) {
                replay_idx = 0;
                continue;
            }
        } else {
            sim->update(Constants::dt / 5, scene);
        }
        this_thread::sleep_for(std::chrono::milliseconds(15));

        // Atomically update renderer_scene
        renderer_scene_updated = true;
        renderer_scene_mutex.lock();
        renderer_scene = scene;
        renderer_scene_mutex.unlock();

        if (save_replay) {
            scene.save(replay_file);
        }

        if (renderer_done) break;
    }

    cout << "Terminating program" << endl;
    glfwTerminate();
    return 0;
}
