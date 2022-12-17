#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include "Constants.hpp"
#include "Scene.hpp"
#include "Simulator.hpp"
#include <iostream>
#include "Renderer.hpp"
#include <memory>
#include "Renderer.hpp"
#include <thread>
#include <chrono>
#include <mutex>
#include <atomic>
#include "ReplayManagers.hpp"

using namespace std;

Scene renderer_scene;
mutex renderer_scene_mutex;
atomic_bool renderer_scene_updated(true);
atomic_bool renderer_done(false);

void render() {
    Renderer renderer;
    static Scene current_scene;

    while (true) {
        // atomically copy renderer_scene
        if (renderer_scene_updated) {
            renderer_scene_mutex.lock();
            current_scene = renderer_scene;
            renderer_scene_mutex.unlock();
            renderer_scene_updated = false;
        }

        // render the copied scene
        if (!renderer.draw(current_scene)) break;
    }

    renderer_done = true;
}

int main(int argc, char* argv[]) {
    cout << "Starting main" << endl;
    // todo: parse these arguments
    bool parallel = true;
    string scene_name = "10000_random_small";
    bool benchmark = false;
    int num_iterations = 100000;
    bool save_replay = true;
    bool play_replay = false;

    if (save_replay && play_replay) {
        cout << "Cannot both save and play replay" << endl;
        return -1;
    }

    string scene_file = string(SCENE_DIR) + scene_name + ".txt";
    cout << "Loading scene from " << scene_file << endl;
    Scene scene(scene_file);
    cout << scene.particles.size() << " particles loaded" << endl;

    cout << "Initiating sim" << endl;
    unique_ptr<Simulator> sim;
    if (parallel) sim.reset(new ParallelSimulator(scene));
    else sim.reset(new SequentialSimulator(scene));

    cout << "Initiating replay managers" << endl;
    string replay_file = string(REPLAY_DIR) + scene_name + ".replay";
    unique_ptr<ReplayWriter> rw;
    if (save_replay) rw.reset(new ReplayWriter(replay_file));
    unique_ptr<ReplayReader> rr;
    if (play_replay) rr.reset(new ReplayReader(replay_file));
    
    cout << "Starting renderer" << endl;
    renderer_scene = scene;
    thread render_thread(render);

    cout << "Starting main loop" << endl;
    for (int i = 0; i < num_iterations; i++) {
        /* Render here */
        if (play_replay) {
            int err = rr->read_scene(&scene);
            if (err) {
                rr->reset();
                continue;
            }
            this_thread::sleep_for(std::chrono::milliseconds(15));
        } else {
            auto st = chrono::high_resolution_clock::now().time_since_epoch();
            sim->update(Constants::dt, scene);
            auto en = chrono::high_resolution_clock::now().time_since_epoch();
            cout << "Update took " << chrono::duration_cast<chrono::milliseconds>(en - st).count() << " ms" << endl;
        }

        // Atomically update renderer_scene
        renderer_scene_mutex.lock();
        renderer_scene = scene;
        renderer_scene_mutex.unlock();
        renderer_scene_updated = true;

        if (save_replay) {
            rw->write_scene(scene);
        }

        if (renderer_done) break;
    }

    cout << "Terminating program" << endl;
    glfwTerminate();
    return 0;
}
